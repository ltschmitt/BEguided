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
beguided.buchholzlab.org {
    # Static assets (matched by path) are the only responses Caddy compresses:
    # the app gzips its own large /results JSON, and SSE (/logs/stream) must
    # never be compressed or the live log stream stalls.
    @static path /static/*
    encode @static gzip zstd

    reverse_proxy 127.0.0.1:8000 {
        header_up X-Real-IP {remote_host}
        header_up X-Forwarded-For {remote_host}
        header_up X-Forwarded-Proto {scheme}
        header_up Host {host}
    }

    # Long-lived caching for the static assets served under /static.
    header @static Cache-Control "public, max-age=86400"

    # Size-rotated access log; no external logrotate needed.
    log {
        output file /var/log/caddy/beguided-access.log {
            roll_size 10MiB
            roll_keep 5
            roll_keep_for 720h
        }
    }
}
```

Caddy's automatic Let's Encrypt issuance kicks in once DNS resolves. If you serve the app under a path prefix (e.g. `example.org/beguided`) instead, wrap the `reverse_proxy` in `handle_path /beguided* { … }` — `handle_path` strips the prefix so FastAPI still sees paths at root. The app itself gzips the large `/results` and `/download/results` payloads (see `SelectiveGZipMiddleware`), so compression works even without Caddy.

### Capacity and parallelism

| Variable | Default (compose) | Meaning |
|---|---|---|
| `BEGUIDED_WORKERS` | `1` | Parallel pipeline worker threads. Threads share the preloaded genome/GTF, so each extra worker costs CPU time but **not** RAM. Set to `cores - 1` (e.g. `3` on a 4-core box). |
| `BEGUIDED_QUEUE_MAX` | `16` | Max queued jobs; further submits get HTTP 503. |
| `BEGUIDED_COOLDOWN_SECONDS` | `15` | Per-IP delay between submits (HTTP 429). |
| `BEGUIDED_MAX_UPLOAD_MB` | `1` | Request body cap (target tables are tiny text). |
| `BEGUIDED_MAX_INPUT_TARGETS` | `500` | Max target rows per upload; also shown in the UI (read from `/config`). |

Do **not** raise `BEGUIDED_WORKERS` above the core count, and never pass a worker count to uvicorn — that spawns separate processes that each reload the genome into RAM.

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
| **Usage statistics** | Anonymous aggregate CSV (job counts, type, target counts); no IP, no user-agent, date granularity only — outside GDPR scope, retained indefinitely |

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
| Usage statistics CSV | Anonymous; retained indefinitely (outlives the 24h job purge) |
| Uploaded custom base editor profiles | Server: in-memory only, lost on restart. Browser: cached in `localStorage` (`beguidedCustomProfiles.v2`) and re-POSTed on every page load so the same browser rehydrates the server copy automatically. Switching browsers or computers starts fresh; clearing site data wipes the local copy. Uploads accept XLSX (recommended, matches the downloaded template) and TSV; the server normalises everything to TSV in memory so the pipeline never sees XLSX bytes. |

The retention sweeper runs every 15 min by default.

## Usage statistics

The web server can record a small, anonymous row for every finished job so
the operator can answer long-term questions — is the service used at all,
how often, and do people run one or two targets or full screens — without
collecting anything that identifies a visitor. Only **web** jobs are
recorded; the CLI pipeline is never counted.

Recording is **on by default**. It never affects job processing: any write
failure is logged and ignored.

### What is recorded

One row per finished job in a monthly CSV (`<stem>_YYYY-MM.csv`):

| Column | Meaning |
|---|---|
| `date` | Server-local calendar day, `YYYY-MM-DD` |
| `mode` | `spreadsheet` or `sequence` |
| `n_targets` | Number of target rows / sequences submitted |
| `n_base_editors` | How many base editors were searched with (just the count) |
| `duration_s` | Wall-clock run time, one decimal |
| `status` | `complete`, `error`, or `cancelled` |

### What is deliberately *not* recorded

No IP address, no user-agent, no job/session identifier, no sub-day
timestamp, and no base-editor **names** (custom profile names are
user-typed free text and are never written). The recorded data therefore
contains no personal data under the GDPR (Recital 26): no legal basis, no
consent, no cookie banner, and indefinite retention are all fine. **The
one rule for future changes: never add an IP address, user-agent, or
sub-day timestamp to this file, or it becomes personal data and needs DPO
review.**

### Configuration

| Variable | Default | Purpose |
|---|---|---|
| `BEGUIDED_USAGE_STATS_PATH` | `/var/lib/beguided/usage_stats.csv` | Base CSV path; written monthly as `<stem>_YYYY-MM.csv` |
| `BEGUIDED_USAGE_STATS_ENABLED` | `true` | Set to `false` to disable recording |

CLI equivalents: `--usage-stats-path` and `--no-usage-stats`.

### Getting the CSV out of Docker

The default path lives on the `beguided-jobs` volume (the same one as
`jobs.db`), so it **survives `docker compose build` + `up -d`** and is
only lost on `docker compose down -v`.

```bash
# Copy the monthly files to the host:
docker compose cp web:/var/lib/beguided/usage_stats_2026-10.csv .
# Or read them in place:
docker compose exec web cat /var/lib/beguided/usage_stats_2026-10.csv
```

To keep the CSVs directly on the host, uncomment the `./usage_stats`
bind mount in `docker-compose.yml` and set
`BEGUIDED_USAGE_STATS_PATH=/var/lib/beguided/usage/usage_stats.csv`.

### Summarizing

```bash
python scripts/summarize_usage.py
python scripts/summarize_usage.py --glob '/var/lib/beguided/usage_stats_*.csv' --last 6
```

Prints per-month job counts, status/mode splits, target-count buckets
(`1-2 / 3-10 / 11-50 / 51-100 / >100`), duration stats, and the average
number of base editors per job.

## Pre-launch checklist

- [ ] `configs/web-legal.yaml` filled in by someone who can read German legalese; reviewed for accuracy
- [ ] `BEGUIDED_ADMIN_TOKEN` and `BEGUIDED_LOG_IP_SALT` set to 32+ random bytes
- [ ] `BEGUIDED_TRUSTED_PROXY_IPS=127.0.0.1,::1` (or Caddy's container IP)
- [ ] Caddy reaches Let's Encrypt and serves HTTPS
- [ ] Browser dev-tools: no requests to non-`buchholzlab.org` hosts on first load; all security headers present
- [ ] `curl -I https://beguided.buchholzlab.org/` returns 200 with the expected headers
- [ ] `curl -I -H "Host: evil.com" https://beguided.buchholzlab.org/` returns 400
- [ ] `pixi run test` is green
- [ ] gitleaks secret scan is green (`pre-commit run --all-files` or CI)
- [ ] `/etc/docker/daemon.json` on the host caps BuildKit cache (see [Docker disk management](#docker-disk-management))
