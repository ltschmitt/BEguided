"use strict";

// DOM elements
const fileInput = document.getElementById("fileInput");
const presetSelect = document.getElementById("presetSelect");
const runBtn = document.getElementById("runBtn");
const cancelBtn = document.getElementById("cancelBtn");
const statusBadge = document.getElementById("statusBadge");
const compatibilityInfo = document.getElementById("compatibilityInfo");
const limitsInfo = document.getElementById("limitsInfo");
const checkAlert = document.getElementById("checkAlert");
const inputAlert = document.getElementById("inputAlert");
const inputPanel = document.getElementById("inputPanel");
const logOutput = document.getElementById("logOutput");
const clearLogsBtn = document.getElementById("clearLogs");
const resultsPanel = document.getElementById("resultsPanel");
const resultsToolbar = document.getElementById("resultsToolbar");
const resultsPlaceholder = document.getElementById("resultsPlaceholder");
const resultsSummary = document.getElementById("resultsSummary");
const bystanderFilter = document.getElementById("bystanderFilter");
const scoreRange = document.getElementById("scoreRange");
const scoreInput = document.getElementById("scoreInput");
const filterSummary = document.getElementById("filterSummary");
const downloadResults = document.getElementById("downloadResults");
const downloadResultsXlsx = document.getElementById("downloadResultsXlsx");
const tabButtons = document.querySelectorAll(".tab-button");
const tabPanels = document.querySelectorAll(".tab-panel");
const toastContainer = document.getElementById("toastContainer");

let uploadedRawB64 = null;
let currentJobId = null;
let currentRevokeToken = null;
let eventSource = null;
let pollHandle = null;
let genomePollHandle = null;
let gtfPollHandle = null;
let currentRecords = [];

// Fetch server config + preset table on load (BEguidedOFF availability, scheduler limits)
(async function init() {
    try {
        const res = await fetch("/config");
        const data = await res.json();
        if (data.ok && data.queue_max != null) {
            // Combined limits summary in the info row.
            limitsInfo.textContent =
                `Queue ${data.queue_max} · ` +
                `${data.cooldown_seconds}s cooldown · ` +
                `${data.max_input_targets} targets · ` +
                `${data.max_upload_mb} MB upload`;
        }
    } catch (err) {
        console.warn("Failed to fetch server config:", err);
    }

    // Genome preload status: one-shot fetch, then poll while loading.
    try {
        await refreshGenomeStatus();
    } catch (err) {
        console.warn("Failed to fetch genome status:", err);
    }
    // GTF preload status: one-shot fetch, then poll while loading.
    try {
        await refreshGtfStatus();
    } catch (err) {
        console.warn("Failed to fetch gtf status:", err);
    }
    // Load the condensed preset table for the Guide tab.
    try {
        const res = await fetch("/preset-table");
        const data = await res.json();
        const container = document.getElementById("presetTableContainer");
        if (data.ok && container) {
            container.innerHTML = renderPresetTable(data.rows);
        } else if (container) {
            container.innerHTML =
                `<p style="color:var(--danger);">` +
                `Failed to load presets: ${escapeHtml(data.error || "unknown error")}` +
                `</p>`;
        }
    } catch (err) {
        const container = document.getElementById("presetTableContainer");
        if (container) {
            container.innerHTML =
                `<p style="color:var(--danger);">Failed to load presets: ${escapeHtml(String(err))}</p>`;
        }
    }
})();

// Toast helper for transient errors (429, 503, network).
function showToast(kind, msg, ttlMs = 8000) {
    const div = document.createElement("div");
    div.className = "alert " + kind;
    div.style.marginBottom = "8px";
    div.style.boxShadow = "0 2px 8px rgba(0,0,0,0.2)";
    div.textContent = msg;
    toastContainer.appendChild(div);
    setTimeout(() => {
        if (div.parentNode) div.parentNode.removeChild(div);
    }, ttlMs);
}

// Genome preload status badge. Polls /genome/status while the preload is in
// flight, then stops on a terminal state (ready / error / disabled).
async function refreshGenomeStatus() {
    const badge = document.getElementById("genomeStatus");
    if (!badge) return;

    let data;
    try {
        const res = await fetch("/genome/status");
        data = await res.json();
    } catch (err) {
        badge.textContent = "genome: status unavailable";
        badge.className = "badge error";
        return;
    }

    if (!data.ok || !data.preload_enabled) {
        badge.textContent = "genome: per-job";
        badge.className = "badge muted";
        if (genomePollHandle) {
            clearInterval(genomePollHandle);
            genomePollHandle = null;
        }
        return;
    }

    if (data.status === "loading") {
        badge.textContent = "genome: loading\u2026";
        badge.className = "badge loading";
        if (!genomePollHandle) {
            genomePollHandle = setInterval(refreshGenomeStatus, 2000);
        }
        return;
    }

    if (data.status === "ready") {
        const label = (data.label || "").trim();
        badge.textContent = label ? `genome: ${label}` : "genome: ready";
        badge.className = "badge ready";
        badge.title = data.path || "";
        if (genomePollHandle) {
            clearInterval(genomePollHandle);
            genomePollHandle = null;
        }
        return;
    }

    if (data.status === "error") {
        badge.textContent = "genome: load failed";
        badge.className = "badge error";
        badge.title = data.error || "";
        if (genomePollHandle) {
            clearInterval(genomePollHandle);
            genomePollHandle = null;
        }
        return;
    }

    // not_loaded or anything unexpected — keep quiet but don't poll.
    badge.textContent = "genome: idle";
    badge.className = "badge muted";
}

// GTF preload status badge. Polls /gtf/status while the preload is in
// flight, then stops on a terminal state (ready / error / disabled).
async function refreshGtfStatus() {
    const badge = document.getElementById("gtfStatus");
    if (!badge) return;

    let data;
    try {
        const res = await fetch("/gtf/status");
        data = await res.json();
    } catch (err) {
        badge.textContent = "gtf: status unavailable";
        badge.className = "badge error";
        return;
    }

    if (!data.ok || !data.preload_enabled) {
        badge.textContent = "gtf: per-job";
        badge.className = "badge muted";
        if (gtfPollHandle) {
            clearInterval(gtfPollHandle);
            gtfPollHandle = null;
        }
        return;
    }

    if (data.status === "loading") {
        badge.textContent = "gtf: loading\u2026";
        badge.className = "badge loading";
        if (!gtfPollHandle) {
            gtfPollHandle = setInterval(refreshGtfStatus, 2000);
        }
        return;
    }

    if (data.status === "ready") {
        const chromCount = data.chromosome_count || 0;
        badge.textContent = `gtf: ready (${chromCount} chrom)`;
        badge.className = "badge ready";
        badge.title = data.path || "";
        if (gtfPollHandle) {
            clearInterval(gtfPollHandle);
            gtfPollHandle = null;
        }
        return;
    }

    if (data.status === "error") {
        badge.textContent = "gtf: load failed";
        badge.className = "badge error";
        badge.title = data.error || "";
        if (gtfPollHandle) {
            clearInterval(gtfPollHandle);
            gtfPollHandle = null;
        }
        return;
    }

    // not_loaded or anything unexpected — keep quiet but don't poll.
    badge.textContent = "gtf: idle";
    badge.className = "badge muted";
}

// Tab switching
tabButtons.forEach(btn => {
    btn.addEventListener("click", () => switchTab(btn.dataset.tab));
});

function switchTab(tabId) {
    tabButtons.forEach(b => {
        b.classList.toggle("active", b.dataset.tab === tabId);
    });
    tabPanels.forEach(p => {
        p.classList.toggle("active", p.id === tabId);
    });
}

async function runCompatibilityCheck() {
    if (!uploadedRawB64) return;
    const form = new FormData();
    form.append("raw_b64", uploadedRawB64);
    form.append("preset", presetSelect.value);
    try {
        const res = await fetch("/check", { method: "POST", body: form });
        const data = await res.json();
        if (!data.ok) return;

        if (data.can_run) {
            checkAlert.style.display = "none";
            compatibilityInfo.textContent = `${data.matched}/${data.total} targets match preset(s)`;
            compatibilityInfo.style.color = "var(--success)";
            runBtn.disabled = false;
        } else {
            showAlert(checkAlert, "error",
                `${data.mismatched}/${data.total} target(s) do not match the selected presets. ` +
                `Supported types: ${data.preset_types.join(", ")}. ` +
                `Please choose a different preset or fix your input.`
            );
            compatibilityInfo.textContent = `${data.matched}/${data.total} targets match preset(s)`;
            compatibilityInfo.style.color = "var(--danger)";
            runBtn.disabled = true;
        }
    } catch (err) {
        showToast("error", `Compatibility check failed: ${err}`);
    }
}

// File upload
fileInput.addEventListener("change", async () => {
    const file = fileInput.files[0];
    if (!file) return;

    const form = new FormData();
    form.append("file", file);

    let res, data;
    try {
        res = await fetch("/upload", { method: "POST", body: form });
        data = await res.json();
    } catch (err) {
        showToast("error", `Upload failed: ${err}`);
        return;
    }

    if (res.status === 413) {
        showToast("error", data.error || "File too large.");
        return;
    }

    if (!data.ok) {
        showAlert(inputAlert, "error", data.error);
        uploadedRawB64 = null;
        runBtn.disabled = true;
        compatibilityInfo.textContent = "";
        checkAlert.style.display = "none";
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${data.error}</p>`;
        return;
    }

    showAlert(inputAlert, "success", `Loaded ${data.row_count} rows. Showing up to 100.`);
    uploadedRawB64 = data.raw_b64;
    inputPanel.innerHTML = renderTable(data.preview, data.columns);
    await runCompatibilityCheck();
    // Switch to the Input tab so the user sees their data.
    switchTab("inputTab");
});

// Re-check when preset changes
presetSelect.addEventListener("change", async () => {
    if (uploadedRawB64) {
        await runCompatibilityCheck();
    }
});

// Run pipeline
runBtn.addEventListener("click", async () => {
    if (!uploadedRawB64) return;

    runBtn.disabled = true;
    setStatus("submitting");
    clearLogs();
    currentRecords = [];
    resultsPanel.querySelectorAll(".results-table-wrap").forEach(n => n.remove());
    if (resultsPlaceholder) {
        resultsPlaceholder.style.display = "block";
        resultsPlaceholder.textContent = "Waiting for pipeline to complete...";
    }
    resultsToolbar.style.display = "none";
    if (resultsSummary) {
        resultsSummary.innerHTML = "";
        resultsSummary.style.display = "none";
    }
    filterSummary.textContent = "Showing 0 of 0";
    downloadResults.style.display = "none";
    downloadResultsXlsx.style.display = "none";
    cancelBtn.style.display = "none";

    // Switch to Logging tab
    tabButtons.forEach(b => b.classList.remove("active"));
    tabPanels.forEach(p => p.classList.remove("active"));
    document.querySelector('[data-tab="logTab"]').classList.add("active");
    document.getElementById("logTab").classList.add("active");

    const form = new FormData();
    form.append("raw_b64", uploadedRawB64);
    form.append("preset", presetSelect.value);

    let res, data;
    try {
        res = await fetch("/run", { method: "POST", body: form });
        data = await res.json();
    } catch (err) {
        setStatus("error");
        logLine("Error starting pipeline: " + err);
        runBtn.disabled = false;
        return;
    }

    if (res.status === 429) {
        showToast("warning",
            `${data.error || "Please wait before submitting again."} ` +
            `(Retry-After: ${res.headers.get("Retry-After") || data.retry_after || "?"}s)`
        );
        setStatus("idle");
        runBtn.disabled = false;
        return;
    }
    if (res.status === 503) {
        showToast("error",
            data.error || "Server busy. Try again in a few seconds."
        );
        setStatus("idle");
        runBtn.disabled = false;
        return;
    }
    if (!data.ok) {
        setStatus("error");
        logLine("Error starting pipeline: " + (data.error || "unknown"));
        runBtn.disabled = false;
        return;
    }

    currentJobId = data.job_id;
    currentRevokeToken = data.revoke_token;
    cancelBtn.style.display = "inline-flex";
    setStatus("queued");
    startLogStream(currentJobId);
    pollStatus(currentJobId);
});

// Cancel button
cancelBtn.addEventListener("click", async () => {
    if (!currentJobId || !currentRevokeToken) return;
    cancelBtn.disabled = true;
    const form = new FormData();
    form.append("job_id", currentJobId);
    form.append("revoke_token", currentRevokeToken);
    try {
        const res = await fetch("/cancel", { method: "POST", body: form });
        const data = await res.json();
        if (!data.ok) {
            showToast("error", `Cancel failed: ${data.error}`);
        } else {
            logLine(`Cancel requested for job ${currentJobId} (status=${data.status}).`);
        }
    } catch (err) {
        showToast("error", `Cancel request failed: ${err}`);
    } finally {
        cancelBtn.disabled = false;
    }
});

// Clear logs
clearLogsBtn.addEventListener("click", clearLogs);

function clearLogs() {
    logOutput.textContent = "";
}

function logLine(text) {
    logOutput.textContent += text + "\n";
    logOutput.scrollTop = logOutput.scrollHeight;
}

function setStatus(status) {
    statusBadge.textContent = status;
    statusBadge.className = "badge " + status;
}

function showAlert(el, kind, msg) {
    el.className = "alert " + kind;
    el.textContent = msg;
    el.style.display = "block";
}

function renderTable(rows, columns, displayColumns = null, groupKeyFn = null) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);padding:16px;">No data to display.</p>`;
    }
    const headers = displayColumns || columns;
    let html = `<table class="data-table"><thead><tr>`;
    headers.forEach(c => { html += `<th>${escapeHtml(c)}</th>`; });
    html += `</tr></thead><tbody>`;
    let prevKey = null;
    let stripeToggle = false;
    rows.forEach(row => {
        let cls = "";
        if (groupKeyFn) {
            const key = groupKeyFn(row);
            if (key !== prevKey) {
                stripeToggle = !stripeToggle;
                prevKey = key;
            }
            cls = stripeToggle ? "target-a" : "target-b";
        }
        html += cls ? `<tr class="${cls}">` : `<tr>`;
        columns.forEach(c => {
            const v = row[c] !== undefined ? String(row[c]) : "";
            html += `<td>${escapeHtml(v)}</td>`;
        });
        html += `</tr>`;
    });
    html += `</tbody></table>`;
    return html;
}

function escapeHtml(str) {
    const div = document.createElement("div");
    div.textContent = str;
    return div.innerHTML;
}

function displayColumnName(name) {
    return name.replace(/_/g, " ");
}

function renderPresetTable(rows) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);">No presets available.</p>`;
    }
    const cols = ["name", "editing", "general", "high", "pam", "cas"];
    const labels = {
        name: "Name",
        editing: "Editing",
        general: "General window",
        high: "High-activity window",
        pam: "PAM",
        cas: "Cas",
    };
    let html = `<table class="data-table"><thead><tr>`;
    cols.forEach(c => {
        html += `<th>${escapeHtml(labels[c])}</th>`;
    });
    html += `</tr></thead><tbody>`;
    rows.forEach(row => {
        html += `<tr>`;
        cols.forEach(c => {
            const v = row[c] !== undefined ? String(row[c]) : "";
            html += `<td>${escapeHtml(v)}</td>`;
        });
        html += `</tr>`;
    });
    html += `</tbody></table>`;
    return html;
}

// SSE log streaming
function startLogStream(jobId) {
    if (eventSource) { eventSource.close(); }

    eventSource = new EventSource("/logs/stream?job_id=" + encodeURIComponent(jobId));
    eventSource.onmessage = (ev) => {
        logLine(ev.data);
    };
    eventSource.onerror = () => {
        // Auto-reconnect or let server close naturally when done
    };
}

// Status polling
async function pollStatus(jobId) {
    const check = async () => {
        try {
            const res = await fetch("/status?job_id=" + encodeURIComponent(jobId));
            const data = await res.json();

            // Update badge for queued/running transitions.
            if (data.status === "running" || data.status === "queued" || data.status === "cancelled") {
                setStatus(data.status);
            }

            if (data.status === "complete" || data.status === "error" || data.status === "cancelled") {
                setStatus(data.status);
                if (eventSource) { eventSource.close(); eventSource = null; }
                runBtn.disabled = false;
                cancelBtn.style.display = "none";

                if (data.status === "complete") {
                    await fetchResults(jobId);
                } else if (data.status === "cancelled") {
                    showToast("info", "Run was cancelled.");
                } else {
                    showToast("error", "Pipeline failed. Check logs for details.");
                }
                return;
            }

            pollHandle = setTimeout(check, 1500);
        } catch (err) {
            // transient network error; try again
            pollHandle = setTimeout(check, 3000);
        }
    };
    pollHandle = setTimeout(check, 1500);
}

// Fetch results
async function fetchResults(jobId) {
    const res = await fetch("/results?job_id=" + encodeURIComponent(jobId));
    const data = await res.json();

    if (!data.ok || !data.records || data.records.length === 0) {
        currentRecords = [];
        if (resultsPlaceholder) {
            resultsPlaceholder.style.display = "block";
            resultsPlaceholder.textContent = data.message || "No results returned.";
        }
        resultsToolbar.style.display = "none";
        if (resultsSummary) {
            resultsSummary.innerHTML = "";
            resultsSummary.style.display = "none";
        }
        showToast("info", data.message || "No results returned.");
        return;
    }

    currentRecords = data.records;
    if (resultsPlaceholder) resultsPlaceholder.style.display = "none";
    resultsToolbar.style.display = "flex";
    renderResultsTable();

    // Auto-switch to results tab
    tabButtons.forEach(b => b.classList.remove("active"));
    tabPanels.forEach(p => p.classList.remove("active"));
    document.querySelector('[data-tab="resultsTab"]').classList.add("active");
    document.getElementById("resultsTab").classList.add("active");
}

// Dynamic filtering for the Results tab. Mirrors beguided.filtering.filters.
function getActiveFilters() {
    let score = parseFloat(scoreInput.value);
    if (!Number.isFinite(score)) score = 0.0;
    if (score < 0) score = 0;
    if (score > 1) score = 1;
    return {
        mode: bystanderFilter.value,
        score,
    };
}

function filterRecords(records, mode, minScore) {
    return records.filter(row => {
        // On-target score filter. Missing scores are kept.
        if (minScore > 0) {
            const raw = row.guide_score;
            if (raw !== undefined && raw !== null && raw !== "") {
                const num = Number(raw);
                if (Number.isFinite(num) && num < minScore) return false;
            }
        }

        if (mode === "all") return true;

        const effect = (row.bystander_effect ?? "").toString();
        if (mode === "no_bystanders") {
            if (effect) return effect === "no_bystanders";
            const positions = (row.bystander_positions ?? "").toString().trim();
            return positions === "";
        }
        if (mode === "silent") {
            return effect === "silent" || effect === "no_bystanders";
        }
        return true;
    });
}

function renderResultsTable() {
    const wrap = resultsPanel.querySelector(".results-table-wrap");
    if (!wrap) {
        const div = document.createElement("div");
        div.className = "results-table-wrap";
        resultsPanel.appendChild(div);
    }
    const target = resultsPanel.querySelector(".results-table-wrap");

    const { mode, score } = getActiveFilters();
    const filtered = filterRecords(currentRecords, mode, score);

    if (filtered.length === 0) {
        target.innerHTML = `<p style="color:var(--muted);padding:16px;">No guides match the current filters.</p>`;
    } else {
        const columns = Object.keys(filtered[0]);
        const displayColumns = columns.map(displayColumnName);
        target.innerHTML = renderTable(
            filtered,
            columns,
            displayColumns,
            r => `${r.chromosome ?? ""}__${r.position ?? ""}`
        );
    }

    renderResultsSummary(filtered);

    filterSummary.textContent = `Showing ${filtered.length} of ${currentRecords.length}`;

    if (currentJobId && currentRecords.length > 0) {
        const params = new URLSearchParams({
            job_id: currentJobId,
            filter_mode: mode,
            guide_score_threshold: score.toString(),
        });
        const baseHref = "/download/results?" + params.toString();
        downloadResults.href = baseHref;
        downloadResults.style.display = "inline-flex";
        downloadResultsXlsx.href = baseHref + "&format=xlsx";
        downloadResultsXlsx.style.display = "inline-flex";
    }
}

// Pivoted per-preset summary of distinct targets visible under the current
// filters. Two rows: the top header row holds the be_preset names with
// "all" first, the second row holds the corresponding count of distinct
// (chromosome, position) targets.
function renderResultsSummary(filtered) {
    if (!resultsSummary) return;
    if (!filtered || filtered.length === 0) {
        resultsSummary.innerHTML = "";
        resultsSummary.style.display = "none";
        return;
    }

    const perPreset = new Map();
    const allTargets = new Set();
    for (const row of filtered) {
        const preset = (row.be_preset ?? "").toString() || "(unknown)";
        const key = `${row.chromosome ?? ""}__${row.position ?? ""}`;
        if (!perPreset.has(preset)) perPreset.set(preset, new Set());
        perPreset.get(preset).add(key);
        allTargets.add(key);
    }

    const presets = [...perPreset.keys()].sort((a, b) => {
        if (a === "(unknown)") return 1;
        if (b === "(unknown)") return -1;
        return a.localeCompare(b);
    });

    let html = `<table class="data-table summary-table"><thead><tr>`;
    html += `<th>${escapeHtml(displayColumnName("be_preset"))}</th><th>all</th>`;
    for (const p of presets) {
        html += `<th>${escapeHtml(p)}</th>`;
    }
    html += `</tr></thead><tbody>`;
    html += `<tr><th scope="row">targets</th>`;
    html += `<td>${allTargets.size}</td>`;
    for (const p of presets) {
        html += `<td>${perPreset.get(p).size}</td>`;
    }
    html += `</tr></tbody></table>`;

    resultsSummary.innerHTML = html;
    resultsSummary.style.display = "block";
}

// Filter event wiring (two-way slider/number sync).
bystanderFilter.addEventListener("change", () => {
    if (currentRecords.length > 0) renderResultsTable();
});
scoreRange.addEventListener("input", () => {
    scoreInput.value = scoreRange.value;
    if (currentRecords.length > 0) renderResultsTable();
});
scoreInput.addEventListener("input", () => {
    let v = parseFloat(scoreInput.value);
    if (!Number.isFinite(v)) return;
    if (v < 0) v = 0;
    if (v > 1) v = 1;
    scoreRange.value = v.toString();
    if (currentRecords.length > 0) renderResultsTable();
});
