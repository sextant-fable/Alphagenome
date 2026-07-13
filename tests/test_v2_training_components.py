from __future__ import annotations

import unittest
from pathlib import Path
import tempfile

import pyBigWig
import torch

from alphagenome_pytorch.heads import targets_scaling
from alphagenome_pytorch.losses import multinomial_loss
from scripts import v2_training_components as components
from scripts import compute_v2_track_means
from scripts import v2_gpu_resources


class V2TrainingComponentsTest(unittest.TestCase):
    def test_gpu_selector_uses_only_idle_gpu_2_then_3(self) -> None:
        resources = {
            "gpus": [
                {"index": 0, "uuid": "u0", "memory_free_mib": 81000, "utilization_percent": 0},
                {"index": 1, "uuid": "u1", "memory_free_mib": 81000, "utilization_percent": 0},
                {"index": 2, "uuid": "u2", "memory_free_mib": 81000, "utilization_percent": 0},
                {"index": 3, "uuid": "u3", "memory_free_mib": 81000, "utilization_percent": 0},
            ],
            "processes": [{"gpu_uuid": "u2", "pid": 1}],
        }
        self.assertEqual(
            v2_gpu_resources.select_available(resources, count=1), [3]
        )
        resources["processes"] = []
        self.assertEqual(
            v2_gpu_resources.select_available(resources, count=2), [2, 3]
        )

    def test_nonzero_mean_counts_only_covered_bases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "signal.bw"
            with pyBigWig.open(str(path), "w") as bigwig:
                bigwig.addHeader([("I", 10), ("II", 10)])
                bigwig.addEntries(
                    ["I", "I", "II"],
                    [0, 5, 0],
                    ends=[5, 10, 10],
                    values=[2.0, 0.0, 4.0],
                )
            totals = compute_v2_track_means.nonzero_totals(path, ("I", "II"))
            self.assertEqual(totals["I"], (10.0, 5))
            self.assertEqual(totals["II"], (40.0, 10))
            self.assertEqual(
                compute_v2_track_means.mean_for_chromosomes(totals, ["I", "II"]),
                50 / 15,
            )

    def test_scale_unscale_matches_formula_and_round_trips(self) -> None:
        target = torch.tensor([[[0.0, 2.0, 20.0], [1.0, 4.0, 40.0]]])
        means = torch.tensor([2.0, 4.0])
        scaled = components.scale_targets_model_space(target, means, resolution=1)
        normalized = target / means.view(1, 2, 1)
        powered = normalized.pow(0.75)
        expected = torch.where(
            powered > 10.0,
            2.0 * torch.sqrt(powered * 10.0) - 10.0,
            powered,
        )
        torch.testing.assert_close(scaled, expected)
        unscaled = components.unscale_predictions_experimental_space(
            scaled, means, resolution=1
        )
        torch.testing.assert_close(unscaled, target, rtol=1e-5, atol=1e-5)

    def test_dual_loss_uses_scaled_targets_at_both_resolutions(self) -> None:
        torch.manual_seed(7)
        target_1 = torch.rand(1, 3, 1024) * 5
        target_128 = target_1.reshape(1, 3, 8, 128).sum(dim=-1)
        means = torch.tensor([1.5, 2.0, 3.0])
        scaled_1 = components.scale_targets_model_space(target_1, means, 1)
        scaled_128 = components.scale_targets_model_space(target_128, means, 128)
        prediction_1 = scaled_1.clone().requires_grad_(True)
        prediction_128 = scaled_128.clone().requires_grad_(True)
        track_mask = torch.ones(1, 3, 1, dtype=torch.bool)
        track_strand = torch.zeros(1, 3, dtype=torch.int8)
        gene_mask = torch.ones(1, 2, 1024, dtype=torch.bool)
        loss, metrics = components.dual_resolution_paper_loss(
            {1: prediction_1, 128: prediction_128},
            {1: target_1, 128: target_128},
            track_means=means,
            track_mask=track_mask,
            track_strand=track_strand,
            gene_mask=gene_mask,
            gene_weight=0.0,
        )
        expected_1 = multinomial_loss(
            y_true=targets_scaling(
                target_1,
                means.unsqueeze(0),
                resolution=1,
                apply_squashing=True,
                channels_last=False,
            ),
            y_pred=prediction_1,
            mask=track_mask,
            multinomial_resolution=128,
            positional_weight=5.0,
            channels_last=False,
        )["loss"]
        expected_128 = multinomial_loss(
            y_true=targets_scaling(
                target_128,
                means.unsqueeze(0),
                resolution=128,
                apply_squashing=True,
                channels_last=False,
            ),
            y_pred=prediction_128,
            mask=track_mask,
            multinomial_resolution=1,
            positional_weight=5.0,
            channels_last=False,
        )["loss"]
        torch.testing.assert_close(loss, expected_1 + expected_128)
        self.assertTrue(torch.isfinite(metrics["loss"]))
        loss.backward()
        self.assertTrue(torch.isfinite(prediction_1.grad).all())
        self.assertTrue(torch.isfinite(prediction_128.grad).all())

    def test_gene_cross_track_loss_is_strand_aware(self) -> None:
        target = torch.ones(1, 3, 8)
        prediction = target.clone()
        prediction[0, 1, :4] = 1000
        mask = torch.ones(1, 3, 1, dtype=torch.bool)
        strands = torch.tensor([[1, -1, 0]], dtype=torch.int8)
        genes = torch.zeros(1, 2, 8, dtype=torch.bool)
        genes[0, 0, :4] = True
        genes[0, 1, 4:] = True
        unchanged = components.gene_cross_track_loss(
            prediction, target, mask, strands, genes
        )
        torch.testing.assert_close(unchanged, torch.tensor(0.0), atol=1e-6, rtol=0)
        prediction[0, 1, 4:] = 1000
        changed = components.gene_cross_track_loss(
            prediction, target, mask, strands, genes
        )
        self.assertGreater(float(changed), 0.1)

    def test_reverse_complement_twice_restores_all_tensors(self) -> None:
        item = {
            "dna_sequence": torch.arange(32).reshape(4, 8).float(),
            "target_1bp": torch.arange(24).reshape(3, 8).float(),
            "target_128bp": torch.arange(6).reshape(3, 2).float(),
            "track_mask": torch.tensor([[True], [False], [True]]),
            "track_strand": torch.tensor([1, -1, 0], dtype=torch.int8),
            "gene_mask": torch.tensor(
                [[True] * 4 + [False] * 4, [False] * 4 + [True] * 4]
            ),
            "core_mask": torch.tensor([True, True, False, False, False, False, True, True]),
        }
        paired = [1, 0, 2]
        twice = components.reverse_complement_item(
            components.reverse_complement_item(item, paired), paired
        )
        for key, value in item.items():
            torch.testing.assert_close(twice[key], value)

    def test_synchronized_crop_recomputes_128bp_sum(self) -> None:
        target = torch.arange(512, dtype=torch.float32).reshape(1, 512)
        item = {
            "dna_sequence": torch.ones(4, 512),
            "target_1bp": target,
            "target_128bp": target.reshape(1, 4, 128).sum(dim=-1),
            "gene_mask": torch.ones(2, 512, dtype=torch.bool),
            "core_mask": torch.ones(512, dtype=torch.bool),
        }
        cropped = components.synchronized_crop(item, start=128, length=256)
        torch.testing.assert_close(cropped["target_1bp"], target[:, 128:384])
        torch.testing.assert_close(
            cropped["target_128bp"], target[:, 128:384].reshape(1, 2, 128).sum(dim=-1)
        )

    def test_worm_embedding_is_a_real_third_trainable_row(self) -> None:
        class Embedder(torch.nn.Module):
            def __init__(self, dimension: int) -> None:
                super().__init__()
                self.num_organisms = 2
                self.organism_embed = torch.nn.Embedding(2, dimension)

        class FakeModel(torch.nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.num_organisms = 2
                self.organism_embed = torch.nn.Embedding(2, 4)
                self.embedder_128bp = Embedder(5)
                self.embedder_1bp = Embedder(6)
                self.embedder_pair = Embedder(7)

        model = FakeModel()
        original = model.organism_embed.weight.detach().clone()
        paths = components.add_c_elegans_organism_embeddings(model)
        self.assertEqual(model.num_organisms, 3)
        self.assertEqual(len(paths), 4)
        torch.testing.assert_close(model.organism_embed.weight[:2], original)
        torch.testing.assert_close(
            model.organism_embed.weight[2], original.mean(dim=0)
        )
        model.organism_embed(torch.tensor([2])).sum().backward()
        self.assertEqual(int(torch.count_nonzero(model.organism_embed.weight.grad[:2])), 0)
        self.assertGreater(int(torch.count_nonzero(model.organism_embed.weight.grad[2])), 0)

    def test_from_scratch_baseline_forward_backward_is_finite(self) -> None:
        torch.manual_seed(3)
        model = components.WormSequenceBaseline(n_tracks=5, hidden_channels=16)
        dna = torch.rand(1, 4, 1024)
        outputs = model(dna)
        self.assertEqual(outputs[1].shape, (1, 5, 1024))
        self.assertEqual(outputs[128].shape, (1, 5, 8))
        loss = outputs[1].mean() + outputs[128].mean()
        loss.backward()
        self.assertTrue(all(torch.isfinite(parameter.grad).all() for parameter in model.parameters()))

    def test_alphagenome_wrapper_transposes_dna_to_nlc(self) -> None:
        class FakeBase(torch.nn.Module):
            def encode(self, dna, organism, resolutions, channels_last):
                self.observed_shape = tuple(dna.shape)
                self.observed_organism = organism.detach().clone()
                self.observed_resolutions = tuple(resolutions)
                self.observed_channels_last = channels_last
                return {
                    "embeddings_1bp": torch.ones(dna.shape[0], 1536, dna.shape[1]),
                    "embeddings_128bp": torch.ones(dna.shape[0], 3072, dna.shape[1] // 128),
                }

        base = FakeBase()
        model = components.AlphaGenomeRnaModel(
            base,
            n_tracks=3,
            track_means=torch.ones(3),
            base_organism_index=2,
            encode_requires_grad=False,
        )
        outputs = model(torch.ones(1, 4, 256))
        self.assertEqual(base.observed_shape, (1, 256, 4))
        self.assertEqual(base.observed_organism.tolist(), [2])
        self.assertEqual(base.observed_resolutions, (1, 128))
        self.assertFalse(base.observed_channels_last)
        self.assertEqual(outputs[1].shape, (1, 3, 256))
        self.assertEqual(outputs[128].shape, (1, 3, 2))

    def test_random_shift_is_seeded_bounded_and_epoch_dependent(self) -> None:
        class FakeDataset(torch.utils.data.Dataset):
            intervals = [
                {
                    "chromosome": "I",
                    "start": "100",
                    "end": "356",
                    "role": "train",
                }
            ]
            fai = {"I": (500, 0, 0, 0)}

            def __len__(self) -> int:
                return 1

            def get_subwindow_item(
                self,
                index: int,
                *,
                shift_bp: int,
                crop_offset_bp: int,
                crop_length_bp: int,
            ):
                target = torch.ones(2, crop_length_bp)
                return {
                    "dna_sequence": torch.ones(4, crop_length_bp),
                    "target_1bp": target,
                    "target_128bp": target.reshape(
                        2, crop_length_bp // 128, 128
                    ).sum(dim=-1),
                    "track_mask": torch.ones(2, 1, dtype=torch.bool),
                    "track_strand": torch.zeros(2, dtype=torch.int8),
                    "gene_mask": torch.ones(2, crop_length_bp, dtype=torch.bool),
                    "core_mask": torch.ones(crop_length_bp, dtype=torch.bool),
                    "shift_bp": shift_bp,
                    "reverse_complemented": False,
                }

        augmented = components.AugmentedV2Dataset(
            FakeDataset(),
            max_shift_bp=1024,
            reverse_complement_probability=0.0,
            seed=17,
        )
        first = augmented[0]["shift_bp"]
        second = augmented[0]["shift_bp"]
        self.assertEqual(first, second)
        self.assertGreaterEqual(first, -100)
        self.assertLessEqual(first, 144)
        augmented.set_epoch(1)
        third = augmented[0]["shift_bp"]
        self.assertGreaterEqual(third, -100)
        self.assertLessEqual(third, 144)
        self.assertNotEqual(first, third)


if __name__ == "__main__":
    unittest.main()
