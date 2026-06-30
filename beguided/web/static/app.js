"use strict";

// DOM elements
const fileInput = document.getElementById("fileInput");
const presetSelect = document.getElementById("presetSelect");
const filterSelect = document.getElementById("filterSelect");
const scoreThreshold = document.getElementById("scoreThreshold");
const offtargetToggle = document.getElementById("offtargetToggle");
const runBtn = document.getElementById("runBtn");
const statusBadge = document.getElementById("statusBadge");
const fileInfo = document.getElementById("fileInfo");
const compatibilityInfo = document.getElementById("compatibilityInfo");
const checkAlert = document.getElementById("checkAlert");
const inputAlert = document.getElementById("inputAlert");
const inputPanel = document.getElementById("inputPanel");
const logOutput = document.getElementById("logOutput");
const clearLogsBtn = document.getElementById("clearLogs");
const resultsAlert = document.getElementById("resultsAlert");
const resultsPanel = document.getElementById("resultsPanel");
const downloadResults = document.getElementById("downloadResults");
const tabButtons = document.querySelectorAll(".tab-button");
const tabPanels = document.querySelectorAll(".tab-panel");

let uploadedRawB64 = null;
let currentJobId = null;
let eventSource = null;

// Fetch server config on load (BEguidedOFF availability, limits)
(async function init() {
    try {
        const res = await fetch("/config");
        const data = await res.json();
        if (data.ok && !data.beguidedoff_available) {
            offtargetToggle.checked = false;
            offtargetToggle.disabled = true;
            showAlert(checkAlert, "warning",
                "Off-target search is unavailable: BEguidedOFF is not installed."
            );
        }
    } catch (err) {
        console.warn("Failed to fetch server config:", err);
    }
})();

// Tab switching
tabButtons.forEach(btn => {
    btn.addEventListener("click", () => {
        tabButtons.forEach(b => b.classList.remove("active"));
        tabPanels.forEach(p => p.classList.remove("active"));
        btn.classList.add("active");
        document.getElementById(btn.dataset.tab).classList.add("active");
    });
});

async function runCompatibilityCheck() {
    if (!uploadedRawB64) return;
    const form = new FormData();
    form.append("raw_b64", uploadedRawB64);
    form.append("preset", presetSelect.value);
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
}

// File upload
fileInput.addEventListener("change", async () => {
    const file = fileInput.files[0];
    if (!file) return;

    const form = new FormData();
    form.append("file", file);

    const res = await fetch("/upload", { method: "POST", body: form });
    const data = await res.json();

    if (!data.ok) {
        showAlert(inputAlert, "error", data.error);
        uploadedRawB64 = null;
        runBtn.disabled = true;
        fileInfo.textContent = "";
        compatibilityInfo.textContent = "";
        checkAlert.style.display = "none";
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${data.error}</p>`;
        return;
    }

    showAlert(inputAlert, "success", `Loaded ${data.row_count} rows. Showing up to 100.`);
    uploadedRawB64 = data.raw_b64;
    fileInfo.textContent = `${file.name} (${data.row_count} rows)`;
    inputPanel.innerHTML = renderTable(data.preview, data.columns);
    await runCompatibilityCheck();
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
    setStatus("running");
    clearLogs();
    resultsPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Waiting for pipeline to complete...</p>`;
    downloadResults.style.display = "none";
    resultsAlert.style.display = "none";

    // Switch to Logging tab
    tabButtons.forEach(b => b.classList.remove("active"));
    tabPanels.forEach(p => p.classList.remove("active"));
    document.querySelector('[data-tab="logTab"]').classList.add("active");
    document.getElementById("logTab").classList.add("active");

    const form = new FormData();
    form.append("raw_b64", uploadedRawB64);
    form.append("preset", presetSelect.value);
    form.append("filter_mode", filterSelect.value);
    form.append("guide_score_threshold", scoreThreshold.value);
    form.append("search_offtargets", offtargetToggle.checked ? "true" : "false");

    const res = await fetch("/run", { method: "POST", body: form });
    const data = await res.json();

    if (!data.ok) {
        setStatus("error");
        logLine("Error starting pipeline: " + (data.error || "unknown"));
        runBtn.disabled = false;
        return;
    }

    currentJobId = data.job_id;
    startLogStream(currentJobId);
    pollStatus(currentJobId);
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

function renderTable(rows, columns) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);padding:16px;">No data to display.</p>`;
    }
    let html = `<table class="data-table"><thead><tr>`;
    columns.forEach(c => { html += `<th>${escapeHtml(c)}</th>`; });
    html += `</tr></thead><tbody>`;
    rows.forEach(row => {
        html += `<tr>`;
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
        const res = await fetch("/status?job_id=" + encodeURIComponent(jobId));
        const data = await res.json();

        if (data.status === "complete" || data.status === "error") {
            setStatus(data.status);
            if (eventSource) { eventSource.close(); eventSource = null; }
            runBtn.disabled = false;

            if (data.status === "complete") {
                await fetchResults(jobId);
            } else {
                showAlert(resultsAlert, "error", "Pipeline failed. Check logs for details.");
            }
            return;
        }

        setTimeout(check, 1500);
    };
    setTimeout(check, 1500);
}

// Fetch results
async function fetchResults(jobId) {
    const res = await fetch("/results?job_id=" + encodeURIComponent(jobId));
    const data = await res.json();

    if (!data.ok || !data.records || data.records.length === 0) {
        resultsPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">No results returned.</p>`;
        showAlert(resultsAlert, "info", data.message || "No results.");
        return;
    }

    const columns = Object.keys(data.records[0]);
    resultsPanel.innerHTML = renderTable(data.records, columns);
    showAlert(resultsAlert, "success", data.message || "Pipeline complete.");
    downloadResults.href = "/download/results?job_id=" + encodeURIComponent(jobId);
    downloadResults.style.display = "inline-flex";

    // Auto-switch to results tab
    tabButtons.forEach(b => b.classList.remove("active"));
    tabPanels.forEach(p => p.classList.remove("active"));
    document.querySelector('[data-tab="resultsTab"]').classList.add("active");
    document.getElementById("resultsTab").classList.add("active");
}
