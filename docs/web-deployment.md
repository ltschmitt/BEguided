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

### Docker disk management

Repeated `docker compose build` invocations can let the BuildKit cache (`/data/docker/overlay2/`) grow without bound — the cache from the previous build is kept until something explicitly evicts it. Cap it once on the host, then forget about it.

Edit `/etc/docker/daemon.json` (create it if absent). If the file already exists, merge the `"buildkit"` key into the existing JSON rather than overwriting:

```json
{
  "features": {
    "buildkit": true
  },
  "buildkit": {
    "gc-keep-storage": "10GB"
  }
}
```

Then restart the daemon:

```bash
sudo systemctl restart docker
```

BuildKit will now garbage-collect its own cache once it crosses the 10 GB ceiling, regardless of how often you rebuild. For an immediate one-shot cleanup:

```bash
docker builder prune -af
```

If the disk is already full from prior builds, run that once after restarting the daemon to reclaim the space.

## GDPR / German-public compliance

Mandatory for any German-public web presence. The audit + remediation history lives in `nogit/compliance_safety.md` (private, gitignored). The operational summary:

| Concern | Mitigation |
|---|---|
| **Impressum** (TMG/DDG §5) | `/impressum` route renders from `configs/web-legal.yaml` |
| **Datenschutzerklärung** (GDPR Art. 13) | `/datenschutz` route renders from the same YAML |
| **Terms of Use** | `/agb` route renders from the same YAML |
| **Cookie banner** | Not needed — the service sets no cookies |
| **Data minimisation** | Pipeline only reads the columns it needs; uploads auto-deleted after 24h; custom profiles persist only in the uploader's browser `localStorage` and the server's RAM |
| **Right to erasure** | `POST /admin/purge-all` with `X-Admin-Token` header purges every job immediately |
| **Special-category data warning** | Footer of every page |
| **Third-party requests** | None. Chart.js is self-hosted with an SRI hash. |
| **Security response headers** | CSP, X-Content-Type-Options, Referrer-Policy, Permissions-Policy, COOP, HSTS (HTTPS only) |
| **Rate limiting** | Per-IP cooldown + queue size cap |
| **Client IP logging** | Hashed with `BEGUIDED_LOG_IP_SALT` (salted SHA-256, first 16 hex chars) |

### Legal YAML

The operator's institutional / DPO / DPA content lives in a gitignored YAML file. The committed template at `configs/web-legal.yaml` is a schema reference with every field blank; it is **not loaded by default**. The operator creates their own file at any path and points the app at it.

The schema is block-based: each legal paragraph (controller, DPO, supervisory authority, accessibility contacts) is **one** list-valued field holding the full multi-line paragraph. Every value may also contain Jinja expressions that reference the other values in the file, so a fact (e.g. the institution name) is defined once and reused everywhere.

```yaml
# /etc/beguided/legal.yaml (example operator file, abridged)
institution_name: "UCC, Institute of X, TU Dresden"
address_lines:
  - "Fetscherstr. 74"
  - "01307 Dresden"
contact_email: "contact@tu-dresden.de"
last_updated: "2026-07-23"

controller_block:
  - "{{ institution_name }}"
  - "{% for line in address_lines %}{{ line }}\n{% endfor %}"
  - "Email: mailto:{{ contact_email }}"

dpo_block:
  - "TU Dresden Data Protection Officer"
  - "Email: mailto:informationssicherheit@tu-dresden.de"

dpa_block:
  - "Saxon Data Protection and Transparency Officer"
  - "Postfach 11 01 32, 01330 Dresden"
```

Block lines are HTML-escaped for output; tokens starting with `mailto:`, `http://` or `https://` are turned into hyperlinks (mailto links display without the `mailto:` prefix). A reference to a variable that is not defined anywhere in the file aborts startup.

Tell the app where to find it (any of):

- env: `BEGUIDED_LEGAL_CONFIG=/etc/beguided/web-legal.yaml`
- CLI flag: `beguided-web --legal-config /etc/beguided/web-legal.yaml`
- docker-compose: `BEGUIDED_LEGAL_CONFIG: "/etc/beguided/web-legal.yaml"` + a bind mount for the file

Behaviour:
- File missing, path empty, or path is a directory → WARNING logged, pages render with empty values. (The directory case covers the docker bind-mount footgun: if the host `legal.yaml` does not exist, Docker creates an empty directory at the mount target inside the container.)
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
| Uploaded custom base editor profiles | Server: in-memory only, lost on restart. Browser: cached in `localStorage` (`beguidedCustomProfiles.v2`) and re-POSTed on every page load so the same browser rehydrates the server copy automatically. Switching browsers or computers starts fresh; clearing site data wipes the local copy. Uploads accept XLSX (recommended, matches the downloaded template) and TSV; the server normalises everything to TSV in memory so the pipeline never sees XLSX bytes. |

The retention sweeper runs every 15 min by default.

## Pre-launch checklist

- [ ] `configs/web-legal.yaml` filled in by someone who can read German legalese; reviewed for accuracy
- [ ] `BEGUIDED_ADMIN_TOKEN` and `BEGUIDED_LOG_IP_SALT` set to 32+ random bytes
- [ ] `BEGUIDED_TRUSTED_PROXY_IPS=127.0.0.1,::1` (or Caddy's container IP)
- [ ] Caddy reaches Let's Encrypt and serves HTTPS
- [ ] Browser dev-tools: no requests to non-`buchholz.org` hosts on first load; all security headers present
- [ ] `curl -I https://www.buchholz.org/beguided/` returns 200 with the expected headers
- [ ] `curl -I -H "Host: evil.com" https://www.buchholz.org/beguided/` returns 400
- [ ] `pixi run test` is green
- [ ] `/etc/docker/daemon.json` on the host caps BuildKit cache (see [Docker disk management](#docker-disk-management))
