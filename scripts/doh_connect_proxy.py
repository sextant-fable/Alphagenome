#!/usr/bin/env python3
"""Localhost-only CONNECT proxy with DNS-over-HTTPS name resolution."""

from __future__ import annotations

from contextlib import contextmanager
import http.client
import json
import select
import socket
import socketserver
import ssl
import threading
import time
from typing import Iterator
from urllib.parse import quote


DOH_ADDRESS = "8.8.8.8"
DOH_HOST = "dns.google"


class DohResolver:
    def __init__(self) -> None:
        self._cache: dict[str, tuple[float, list[str]]] = {}
        self._lock = threading.Lock()

    def resolve(self, hostname: str) -> list[str]:
        try:
            socket.inet_pton(socket.AF_INET, hostname)
            return [hostname]
        except OSError:
            pass
        with self._lock:
            cached = self._cache.get(hostname)
            if cached and cached[0] > time.monotonic():
                return cached[1]

        raw_socket = socket.create_connection((DOH_ADDRESS, 443), timeout=15)
        context = ssl.create_default_context()
        with context.wrap_socket(raw_socket, server_hostname=DOH_HOST) as connection:
            path = f"/resolve?name={quote(hostname, safe='')}&type=A"
            request = (
                f"GET {path} HTTP/1.1\r\n"
                f"Host: {DOH_HOST}\r\n"
                "Accept: application/dns-json\r\n"
                "Connection: close\r\n\r\n"
            )
            connection.sendall(request.encode("ascii"))
            response = http.client.HTTPResponse(connection)
            response.begin()
            if response.status != 200:
                raise RuntimeError(
                    f"DoH lookup failed for {hostname}: HTTP {response.status}"
                )
            payload = json.loads(response.read())
        addresses = [
            answer["data"]
            for answer in payload.get("Answer", [])
            if answer.get("type") == 1
        ]
        if not addresses:
            raise RuntimeError(f"DoH returned no IPv4 address for {hostname}")
        ttl = min(
            (int(answer.get("TTL", 60)) for answer in payload.get("Answer", [])),
            default=60,
        )
        with self._lock:
            self._cache[hostname] = (
                time.monotonic() + max(min(ttl, 300), 15),
                addresses,
            )
        return addresses


class ConnectProxyHandler(socketserver.BaseRequestHandler):
    resolver: DohResolver

    def handle(self) -> None:
        self.request.settimeout(30)
        request = b""
        while b"\r\n\r\n" not in request and len(request) < 65_536:
            chunk = self.request.recv(4096)
            if not chunk:
                return
            request += chunk
        first_line = request.split(b"\r\n", 1)[0].decode("ascii", errors="replace")
        fields = first_line.split()
        if len(fields) != 3 or fields[0].upper() != "CONNECT":
            self.request.sendall(b"HTTP/1.1 405 Method Not Allowed\r\n\r\n")
            return
        target = fields[1]
        hostname, separator, port_text = target.rpartition(":")
        if not separator or not hostname:
            self.request.sendall(b"HTTP/1.1 400 Bad Request\r\n\r\n")
            return
        upstream = None
        try:
            for address in self.resolver.resolve(hostname):
                try:
                    upstream = socket.create_connection(
                        (address, int(port_text)), timeout=30
                    )
                    break
                except OSError:
                    continue
            if upstream is None:
                raise RuntimeError(f"Unable to connect to {target}")
            self.request.sendall(b"HTTP/1.1 200 Connection Established\r\n\r\n")
            self.request.settimeout(None)
            upstream.settimeout(None)
            sockets = [self.request, upstream]
            while True:
                readable, _, _ = select.select(sockets, [], [], 120)
                if not readable:
                    return
                for source in readable:
                    destination = upstream if source is self.request else self.request
                    data = source.recv(65_536)
                    if not data:
                        return
                    destination.sendall(data)
        except Exception:
            try:
                self.request.sendall(b"HTTP/1.1 502 Bad Gateway\r\n\r\n")
            except OSError:
                pass
        finally:
            if upstream is not None:
                upstream.close()


class ThreadingConnectProxy(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


@contextmanager
def local_doh_proxy() -> Iterator[str]:
    resolver = DohResolver()
    handler = type(
        "BoundConnectProxyHandler",
        (ConnectProxyHandler,),
        {"resolver": resolver},
    )
    server = ThreadingConnectProxy(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        host, port = server.server_address
        yield f"http://{host}:{port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
