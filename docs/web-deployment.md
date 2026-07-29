# Web deployment

The web server (`beguided-web`) is a FastAPI app designed to run behind a reverse proxy (Caddy, nginx, Traefik, …). This page covers local use, public deployment, GDPR compliance, and the legal YAML.

## Local use

```bash
# After `pip install -e ".[web]"` or `pixi install --features web`
beguided-web
# → http://localhost:8000
```

For local development with auto-reload:

```bash
beguided-web --reload
```

`scripts/run_web.sh` is a thin wrapper that pins scheduler defaults and forwards extra CLI args; useful for personal development.

## Production deployment with Docker Compose

The repo ships a `docker-compose.yml` that builds the CPU web image and runs it on port 8000. Every deployment value is interpolated from a shell env var with an inline default, so `docker compose up -d` works with no environment file at all.

```bash
docker compose build
docker compose up -d
docker compose logs -f web
```

### Public deployment via Caddy

Caddy terminates TLS in front of the container and reverse-proxies to `127.0.0.1:8000`. The minimum Caddy snippet:

```
www.buchholz.org {
    handle_path /beguided* {
        reverse_proxy 127.0.0.1:8000 {
            header_up X-Real-IP {remote_host}
            header_up X-Forwarded-For {remote_host}
            header_up X-Forwarded-Proto {scheme}
            header_up Host {host}
        }
    }
}
```

`handle_path` strips the `/beguided` prefix so FastAPI sees paths at root — the server is unchanged. Caddy's automatic Let's Encrypt issuance kicks in once DNS resolves.

### Reverse-proxy trust

The web app reads `X-Forwarded-For` / `X-Forwarded-Proto` to recover the real client IP and HTTPS scheme, but only when the upstream peer is on the `BEGUIDED_TRUSTED_PROXY_IPS` allow-list (default `127.0.0.1,::1`). Add Caddy's peer IP if it runs in a separate container.

### Host-header allow-list

`BEGUIDED_ALLOWED_HOSTS` (CSV, default `www.buchholz.org,buchholz.org,localhost,127.0.0.1,testserver`) gates Starlette's `TrustedHostMiddleware`. Requests with any other `Host:` value are rejected with HTTP 400. For production, drop the dev values.

## GDPR / German-public compliance

Mandatory for any German-public web presence. The audit + remediation history lives in `nogit/compliance_safety.md` (private, gitignored). The operational summary:

| Concern | Mitigation |
|---|---|
| **Impressum** (TMG/DDG §5) | `/impressum` route renders from `configs/web-legal.yaml` |
| **Datenschutzerklärung** (GDPR Art. 13) | `/datenschutz` route renders from the same YAML |
| **Terms of Use** | `/agb` route renders from the same YAML |
| **Cookie banner** | Not needed — the service sets no cookies |
| **Data minimisation** | Pipeline only reads the columns it needs; uploads auto-deleted after 24h |
| **Right to erasure** | `POST /admin/purge-all` with `X-Admin-Token` header purges every job immediately |
| **Special-category data warning** | Footer of every page |
| **Third-party requests** | None. Chart.js is self-hosted with an SRI hash. |
| **Security response headers** | CSP, X-Content-Type-Options, Referrer-Policy, Permissions-Policy, COOP, HSTS (HTTPS only) |
| **Rate limiting** | Per-IP cooldown + queue size cap |
| **Client IP logging** | Hashed with `BEGUIDED_LOG_IP_SALT` (salted SHA-256, first 16 hex chars) |

### Legal YAML

The operator's institutional / DPO / DPA content lives in a gitignored YAML file. The committed template at `configs/web-legal.yaml` is empty (every field blank); it is **not loaded by default**. The operator creates their own file at any path and points the app at it.

```yaml
# /etc/beguided/web-legal.yaml (example operator file)
institution_name: "..."
address_lines: ["..."]
vertretungsberechtigter: "..."
contact_email: "..."
contact_phone: "..."
hosting_provider: "..."
content_responsible_mstv: "..."
dpo_name: "..."
dpo_email: "..."
dpo_address: ["..."]
competent_dpa: "..."
log_retention_days: 7
last_updated: "2026-07-23"
```

Tell the app where to find it (any of):

- env: `BEGUIDED_LEGAL_CONFIG=/etc/beguided/web-legal.yaml`
- CLI flag: `beguided-web --legal-config /etc/beguided/web-legal.yaml`
- docker-compose: `BEGUIDED_LEGAL_CONFIG: "/etc/beguided/web-legal.yaml"` + a bind mount for the file

Behaviour:
- File missing or path empty → WARNING logged, pages render with empty values
- File present, malformed YAML → fail-fast at startup
- File present, top-level is not a mapping → fail-fast at startup
- File present, valid → loaded once, every legal page renders with the values

### Compliance secrets

| Variable | Purpose | How to set |
|---|---|---|
| `BEGUIDED_ADMIN_TOKEN` | Gates `POST /admin/purge-all` | `python -c "import secrets; print(secrets.token_urlsafe(32))"` |
| `BEGUIDED_LOG_IP_SALT` | Salts the IP hash in logs | same |

### Retention and data minimisation

| Data | Retention |
|---|---|
| Uploaded target TSV | 24h after job completes (`BEGUIDED_JOB_RETENTION_HOURS`) |
| Generated guide results | 24h after job completes |
| Revoke tokens | Hashed SHA-256 in `jobs.db`; raw token shown to user once |
| Scheduler cooldown map | In-memory, TTL-pruned after 4× `cooldown_seconds` |
| Proxy / web logs | Operator-defined (`legal.log_retention_days`); salted IP |

The retention sweeper runs every 15 min by default.

## Pre-launch checklist

- [ ] `configs/web-legal.yaml` filled in by someone who can read German legalese; reviewed for accuracy
- [ ] `BEGUIDED_ADMIN_TOKEN` and `BEGUIDED_LOG_IP_SALT` set to 32+ random bytes
- [ ] `BEGUIDED_ALLOWED_HOSTS=www.buchholz.org,buchholz.org` (no dev hosts)
- [ ] `BEGUIDED_TRUSTED_PROXY_IPS=127.0.0.1,::1` (or Caddy's container IP)
- [ ] Caddy reaches Let's Encrypt and serves HTTPS
- [ ] Browser dev-tools: no requests to non-`buchholz.org` hosts on first load; all security headers present
- [ ] `curl -I https://www.buchholz.org/beguided/` returns 200 with the expected headers
- [ ] `curl -I -H "Host: evil.com" https://www.buchholz.org/beguided/` returns 400
- [ ] `pixi run test` is green
