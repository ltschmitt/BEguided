"use strict";

// DOM elements
const fileInput = document.getElementById("fileInput");
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
const resultsGuideStatus = document.getElementById("resultsGuideStatus");
const bystanderFilter = document.getElementById("bystanderFilter");
const scoreFilter = document.getElementById("scoreFilter");
const filterSummary = document.getElementById("filterSummary");
const downloadResults = document.getElementById("downloadResults");
const downloadResultsXlsx = document.getElementById("downloadResultsXlsx");
const tabButtons = document.querySelectorAll(".tab-button");
const tabPanels = document.querySelectorAll(".tab-panel");
const subtabButtons = document.querySelectorAll(".subtab-button");
const subtabPanels = document.querySelectorAll(".subtab-panel");
const toastContainer = document.getElementById("toastContainer");
const modeSpreadsheetBtn = document.getElementById("modeSpreadsheet");
const modeSequenceBtn = document.getElementById("modeSequence");
const spreadsheetPanel = document.getElementById("spreadsheetPanel");
const sequencePanel = document.getElementById("sequencePanel");
const sequenceInput = document.getElementById("sequenceInput");
const parseSequenceBtn = document.getElementById("parseSequenceBtn");
const loadSequenceExampleBtn = document.getElementById("loadSequenceExampleBtn");
const loadSpreadsheetExampleBtn = document.getElementById("loadSpreadsheetExampleBtn");
const profilesUpload = document.getElementById("profilesUpload");
const profilesAlert = document.getElementById("profilesAlert");
const profilesSource = document.getElementById("profilesSource");
const profilesTableContainer = document.getElementById("profilesTableContainer");

let currentMode = "sequence";
let uploadedRawB64 = null;
let customGenomeB64 = null;
let currentJobId = null;
let currentRevokeToken = null;
let eventSource = null;
let pollHandle = null;
let genomePollHandle = null;
let gtfPollHandle = null;
let currentRecords = [];
let matchedBaseEditorCount = 0;
let customEditorKey = "custom";
let customEditorLabel = "Custom (uploaded)";
let customProfilesLoaded = false;
let lastProfilesTableSource = null;

// Base editor multi-select state. `availableEditors` is the list of editors
// rendered as checkboxes (fetched from /base-editor-table?source=default).
// `selectedEditors` is the Set of editor names the user has checked.
// Custom profiles are merged into `availableEditors` via source="combined",
// so there is no separate "custom mode" — the multi-select stays active
// whenever any profiles (default or custom) are available.
const beMultiselect = document.getElementById("beMultiselect");
const beMultiselectBtn = document.getElementById("beMultiselectBtn");
const beMultiselectPanel = document.getElementById("beMultiselectPanel");
const beMultiselectLabel = document.getElementById("beMultiselectLabel");
const beCheckboxGroups = document.getElementById("beCheckboxGroups");
let availableEditors = [];
let selectedEditors = new Set();
let customProfilesCount = 0;

// Fetch server config + base editor table on load (BEguidedOFF availability, scheduler limits)
(async function init() {
    try {
        const res = await fetch("/config");
        const data = await res.json();
        if (data.ok && data.max_input_targets != null) {
            // Upload limits summary shown inside the spreadsheet input panel.
            limitsInfo.textContent =
                `Upload a table of genomic targets below (max. ${data.max_input_targets} targets ` +
                `and ${data.max_upload_mb} MB)`;
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
    // Hydrate the custom profiles state and render the Base Editors table.
    try {
        await refreshProfilesStatus({ keepDropdownValue: true });
    } catch (err) {
        console.warn("Failed to fetch profiles status:", err);
    }
    await refreshProfilesTable();

    // Load the default base editor profiles into the multi-select.
    try {
        await loadAvailableEditors();
    } catch (err) {
        console.warn("Failed to load available editors:", err);
    }

    wireMultiselectEvents();
})();

// ---------------------------------------------------------------------------
// Base editor multi-select
// ---------------------------------------------------------------------------

// Fetch the combined (default + custom) profiles table and build the
// grouped checkbox panel.
async function loadAvailableEditors() {
    const res = await fetch(
        `/base-editor-table?source=${encodeURIComponent("combined")}`
    );
    const data = await res.json();
    if (!data.ok) {
        beCheckboxGroups.innerHTML =
            `<p style="color:var(--danger);padding:8px;">Failed to load base editor profiles: ` +
            `${escapeHtml(data.error || "unknown error")}</p>`;
        return;
    }
    availableEditors = data.rows.map(r => ({
        name: r.name,
        editing: r.editing,
        cas: r.cas || "",
        pam: r.pam,
        source: r.source,
    }));
    // Default: all editors selected.
    selectedEditors = new Set(availableEditors.map(r => r.name));
    renderCheckboxGroups();
    updateMultiselectLabel();
    updateCustomToggleVisibility();
}

// Show or hide the "All Custom" quick-toggle based on whether any custom
// editors are present.
function updateCustomToggleVisibility() {
    const btn = beMultiselectPanel.querySelector(
        '.be-quick-toggles button[data-toggle="custom"]'
    );
    if (!btn) return;
    const hasCustom = availableEditors.some(ed => ed.source === "custom");
    btn.style.display = hasCustom ? "" : "none";
}

// Render the checkbox groups. Order:
//   Custom (with optional Cas subgroups) -> ABE -> CBE.
function renderCheckboxGroups() {
    if (availableEditors.length === 0) {
        beCheckboxGroups.innerHTML = "";
        return;
    }

    const customEditors = availableEditors.filter(ed => ed.source === "custom");
    const defaultEditors = availableEditors.filter(ed => ed.source !== "custom");

    let html = "";

    // 1. Custom group (only when there are custom editors).
    if (customEditors.length > 0) {
        html += renderOneGroup("Custom", customEditors, {
            dataSource: "custom",
            headerClass: "be-group-header be-group-header-custom",
        });
    }

    // 2. Default groups (ABE / CBE) — grouped by editing -> cas.
    const editingOrder = ["ABE", "CBE"];
    const groups = new Map();
    for (const ed of defaultEditors) {
        if (!groups.has(ed.editing)) groups.set(ed.editing, new Map());
        const casMap = groups.get(ed.editing);
        if (!casMap.has(ed.cas)) casMap.set(ed.cas, []);
        casMap.get(ed.cas).push(ed);
    }
    const editingKeys = [...groups.keys()].sort((a, b) => {
        const ia = editingOrder.indexOf(a);
        const ib = editingOrder.indexOf(b);
        if (ia === -1 && ib === -1) return a.localeCompare(b);
        if (ia === -1) return 1;
        if (ib === -1) return -1;
        return ia - ib;
    });
    for (const editing of editingKeys) {
        html += renderOneGroup(editing, defaultEditors.filter(
            ed => ed.editing === editing
        ));
    }

    beCheckboxGroups.innerHTML = html;

    // Wire up change listeners on the freshly rendered checkboxes.
    beCheckboxGroups.querySelectorAll('input[type="checkbox"]').forEach(cb => {
        cb.addEventListener("change", onCheckboxChange);
    });
}

// A checkbox was toggled by the user.
function onCheckboxChange(ev) {
    const name = ev.target.dataset.name;
    if (!name) return;
    if (ev.target.checked) {
        selectedEditors.add(name);
    } else {
        selectedEditors.delete(name);
    }
    updateMultiselectLabel();
    onSelectionChanged();
}

// Render a single group block. If any editor in the group has a non-empty
// `cas`, sub-group by Cas. Otherwise list directly under the group header
// with no Cas subgroup.
function renderOneGroup(label, editors, opts = {}) {
    const headerClass = opts.headerClass || "be-group-header";
    const dataSource = opts.dataSource || "";
    const dsAttr = dataSource ? ` data-source="${escapeHtml(dataSource)}"` : "";
    let html = `<div class="be-group"${dsAttr}>`;
    html += `<div class="${escapeHtml(headerClass)}">${escapeHtml(label)}</div>`;

    const useSubgroup = editors.some(ed => ed.cas && ed.cas.length > 0);
    if (useSubgroup) {
        const casMap = new Map();
        for (const ed of editors) {
            const key = ed.cas || "";
            if (!casMap.has(key)) casMap.set(key, []);
            casMap.get(key).push(ed);
        }
        const casKeys = [...casMap.keys()].sort((a, b) => {
            // Empty Cas sorts last so it sits below named Cas variants.
            if (!a) return 1;
            if (!b) return -1;
            return a.localeCompare(b);
        });
        for (const cas of casKeys) {
            html += `<div class="be-subgroup" data-cas="${escapeHtml(cas)}">`;
            if (cas) {
                html += `<div class="be-subgroup-header">${escapeHtml(cas)}</div>`;
            }
            const casEditors = casMap.get(cas).sort((a, b) =>
                a.name.localeCompare(b.name)
            );
            for (const ed of casEditors) {
                html += renderCheckboxRow(ed);
            }
            html += `</div>`;
        }
    } else {
        const sorted = [...editors].sort((a, b) => a.name.localeCompare(b.name));
        for (const ed of sorted) {
            html += renderCheckboxRow(ed);
        }
    }

    html += `</div>`;
    return html;
}

function renderCheckboxRow(ed) {
    let html = `<label class="be-checkbox-row" data-name="${escapeHtml(ed.name)}">`;
    html += `<input type="checkbox" data-name="${escapeHtml(ed.name)}"`;
    if (selectedEditors.has(ed.name)) html += " checked";
    html += `>`;
    html += `<span class="be-name">${escapeHtml(ed.name)}</span>`;
    if (ed.pam) {
        html += `<span class="be-pam">${escapeHtml(ed.pam)}</span>`;
    }
    html += `</label>`;
    return html;
}

// Update the dropdown button label based on the current selection.
function updateMultiselectLabel() {
    const total = availableEditors.length;
    const n = selectedEditors.size;
    if (total === 0) {
        beMultiselectLabel.textContent = "Loading...";
        return;
    }
    if (n === 0) {
        beMultiselectLabel.textContent = "None selected";
    } else if (n === total) {
        beMultiselectLabel.textContent = `All (${total})`;
    } else {
        beMultiselectLabel.textContent = `${n} selected`;
    }
}

// Quick toggle: all / abe / cbe / custom / none.
function applyQuickToggle(toggle) {
    if (toggle === "all") {
        selectedEditors = new Set(availableEditors.map(r => r.name));
    } else if (toggle === "abe") {
        selectedEditors = new Set(
            availableEditors
                .filter(r => r.editing === "ABE")
                .map(r => r.name)
        );
    } else if (toggle === "cbe") {
        selectedEditors = new Set(
            availableEditors
                .filter(r => r.editing === "CBE")
                .map(r => r.name)
        );
    } else if (toggle === "custom") {
        selectedEditors = new Set(
            availableEditors
                .filter(r => r.source === "custom")
                .map(r => r.name)
        );
    } else if (toggle === "none") {
        selectedEditors = new Set();
    }
    // Reflect in the rendered checkboxes.
    beCheckboxGroups.querySelectorAll('input[type="checkbox"]').forEach(cb => {
        cb.checked = selectedEditors.has(cb.dataset.name);
    });
    updateMultiselectLabel();
    onSelectionChanged();
}

// Toggle the dropdown panel open/closed.
function toggleMultiselectPanel(force) {
    const isHidden = beMultiselectPanel.hasAttribute("hidden");
    const shouldOpen = force === true || (force !== false && isHidden);
    if (shouldOpen) {
        beMultiselectPanel.removeAttribute("hidden");
        beMultiselectBtn.setAttribute("aria-expanded", "true");
    } else {
        beMultiselectPanel.setAttribute("hidden", "");
        beMultiselectBtn.setAttribute("aria-expanded", "false");
    }
}

function closeMultiselectPanel() {
    toggleMultiselectPanel(false);
}

// Wire up the multi-select DOM events.
function wireMultiselectEvents() {
    beMultiselectBtn.addEventListener("click", (ev) => {
        ev.stopPropagation();
        if (beMultiselect.classList.contains("disabled")) return;
        toggleMultiselectPanel();
    });

    beMultiselectPanel.addEventListener("click", (ev) => ev.stopPropagation());

    beMultiselectPanel.querySelectorAll(".be-quick-toggles button").forEach(btn => {
        btn.addEventListener("click", () => applyQuickToggle(btn.dataset.toggle));
    });

    document.addEventListener("click", (ev) => {
        if (!beMultiselect.contains(ev.target)) closeMultiselectPanel();
    });

    document.addEventListener("keydown", (ev) => {
        if (ev.key === "Escape") closeMultiselectPanel();
    });
}

// Re-run compatibility check + refresh the Base Editors table whenever the
// selection changes.
function onSelectionChanged() {
    lastProfilesTableSource = null; // force table refresh on next ask
    if (uploadedRawB64) {
        runCompatibilityCheck();
    }
}

// Return the value to send to /check and /run: always a comma-separated
// list of the currently checked editor names (default + custom mixed).
function baseEditorSelectionValue() {
    return [...selectedEditors].sort().join(",");
}

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
        const label = (data.label || "").trim();
        badge.textContent = label ? `gtf: ${label}` : "gtf: ready";
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

// Sub-tab switching (nested tabs inside the Results tab)
subtabButtons.forEach(btn => {
    btn.addEventListener("click", () => switchSubTab(btn.dataset.subtab));
});

function switchSubTab(subtabId) {
    subtabButtons.forEach(b => {
        b.classList.toggle("active", b.dataset.subtab === subtabId);
    });
    subtabPanels.forEach(p => {
        p.classList.toggle("active", p.id === subtabId);
    });
}

// Mode switcher: spreadsheet vs sequence input.
function setMode(mode) {
    if (mode !== "spreadsheet" && mode !== "sequence") return;
    currentMode = mode;
    modeSpreadsheetBtn.classList.toggle("active", mode === "spreadsheet");
    modeSequenceBtn.classList.toggle("active", mode === "sequence");
    spreadsheetPanel.style.display = mode === "spreadsheet" ? "" : "none";
    sequencePanel.style.display = mode === "sequence" ? "" : "none";

    // Move the shared Run button into the active input panel's action row so
    // it stays visually anchored next to the input controls. The 'right' class
    // pushes it to the far edge of the row in spreadsheet mode; in sequence
    // mode it must sit immediately after the Parse button, so we strip it.
    const targetRow = mode === "spreadsheet"
        ? spreadsheetPanel.querySelector(".input-actions-row")
        : sequencePanel.querySelector(".input-meta-row");
    if (targetRow) {
        targetRow.appendChild(runBtn);
    }
    runBtn.classList.toggle("right", mode === "spreadsheet");

    // Reset state when switching modes: a previously uploaded file is
    // irrelevant for sequence input and vice-versa.
    uploadedRawB64 = null;
    customGenomeB64 = null;
    matchedBaseEditorCount = 0;
    runBtn.disabled = true;
    compatibilityInfo.textContent = "";
    checkAlert.style.display = "none";
    inputAlert.style.display = "none";
    inputPanel.innerHTML = `<p style="color: var(--muted);padding:16px;">${
        mode === "spreadsheet"
            ? "Upload a targets file to preview input data."
            : "Parse sequences to preview input data."
    }</p>`;
}

modeSpreadsheetBtn.addEventListener("click", () => setMode("spreadsheet"));
modeSequenceBtn.addEventListener("click", () => setMode("sequence"));

async function runCompatibilityCheck() {
    if (!uploadedRawB64) return;
    // Refuse to check with an empty selection.
    if (selectedEditors.size === 0) {
        runBtn.disabled = true;
        compatibilityInfo.textContent = "Select at least one base editor.";
        compatibilityInfo.style.color = "var(--muted)";
        matchedBaseEditorCount = 0;
        return;
    }
    const form = new FormData();
    form.append("raw_b64", uploadedRawB64);
    form.append("base_editor", baseEditorSelectionValue());
    try {
        const res = await fetch("/check", { method: "POST", body: form });
        const data = await res.json();
        if (!data.ok) return;

        if (data.can_run) {
            checkAlert.style.display = "none";
            compatibilityInfo.textContent = `${data.matched}/${data.total} targets match the selected base editor profile(s)`;
            compatibilityInfo.style.color = "var(--success)";
            matchedBaseEditorCount = data.matched;
            runBtn.disabled = false;
        } else {
            showAlert(checkAlert, "error",
                `${data.mismatched}/${data.total} target(s) do not match the selected base editor profile(s). ` +
                `Supported types: ${data.editor_types.join(", ")}. ` +
                `Please choose a different profile or fix your input.`
            );
            compatibilityInfo.textContent = `${data.matched}/${data.total} targets match the selected base editor profile(s)`;
            compatibilityInfo.style.color = "var(--danger)";
            matchedBaseEditorCount = data.matched;
            runBtn.disabled = true;
        }
    } catch (err) {
        showToast("error", `Compatibility check failed: ${err}`);
    }
}

// Sequence-mode parse handler. Posts the textarea content to /sequence/parse
// and, on success, mirrors the spreadsheet-mode state (raw_b64 + preview +
// compatibility check).
async function parseSequences() {
    const text = sequenceInput.value;
    if (!text.trim()) {
        showToast("error", "Paste at least one sequence.");
        runBtn.disabled = true;
        matchedBaseEditorCount = 0;
        compatibilityInfo.textContent = "";
        checkAlert.style.display = "none";
        return;
    }

    // Switch to the Input tab immediately so the user sees feedback — and any
    // error messages — the moment they click Parse.
    switchTab("inputTab");
    inputAlert.style.display = "none";
    checkAlert.style.display = "none";
    compatibilityInfo.textContent = "";
    inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Parsing sequences\u2026</p>`;

    const form = new FormData();
    form.append("text", text);
    let res, data;
    try {
        res = await fetch("/sequence/parse", { method: "POST", body: form });
        data = await res.json();
    } catch (err) {
        showAlert(inputAlert, "error", `Parse failed: ${err}`);
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Parse failed.</p>`;
        return;
    }

    if (res.status === 413) {
        showAlert(inputAlert, "error", data.error || "Input too large.");
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(data.error || "Input too large.")}</p>`;
        return;
    }

    if (!data.ok) {
        showAlert(inputAlert, "error", data.error);
        uploadedRawB64 = null;
        customGenomeB64 = null;
        matchedBaseEditorCount = 0;
        runBtn.disabled = true;
        compatibilityInfo.textContent = "";
        checkAlert.style.display = "none";
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(data.error)}</p>`;
        return;
    }

    showAlert(inputAlert, "success", `Parsed ${data.row_count} sequence(s). Showing up to 100.`);
    uploadedRawB64 = data.raw_b64;
    customGenomeB64 = data.custom_genome_b64;
    inputPanel.innerHTML = renderSequencePreview(data.preview);
    await runCompatibilityCheck();
}

parseSequenceBtn.addEventListener("click", parseSequences);

// Sequence-mode "Load example" button. Populates the textarea with a small
// multi-record FASTA (matches the example shown in the textarea's explanation
// text) and then runs the regular parse flow.
const SEQUENCE_EXAMPLE_FASTA =
    ">my_target_1\n" +
    "ATAGAGTATAGACCC[A>G]CTATCCAGTCTACCTAGGGAATC\n" +
    ">my_target_2\n" +
    "GGGCATAGCCATGTAGCATT[T>C]GCAATTGGCCATGCATTTACC";

if (loadSequenceExampleBtn) {
    loadSequenceExampleBtn.addEventListener("click", () => {
        if (currentMode !== "sequence") setMode("sequence");
        sequenceInput.value = SEQUENCE_EXAMPLE_FASTA;
        parseSequences();
    });
}

// File upload
fileInput.addEventListener("change", async () => {
    if (currentMode !== "spreadsheet") return;
    const file = fileInput.files[0];
    if (!file) return;
    await uploadTargetsFile(file);
    fileInput.value = "";
});

// Reusable upload pipeline for spreadsheet-mode targets. Used by both the
// <input type="file"> change handler and the "Load example" button (which
// synthesizes a File from in-memory text).
async function uploadTargetsFile(file) {
    if (currentMode !== "spreadsheet") return;

    // Switch to the Input tab immediately so the user sees feedback — and any
    // error messages — the moment they pick a file.
    switchTab("inputTab");
    inputAlert.style.display = "none";
    checkAlert.style.display = "none";
    compatibilityInfo.textContent = "";
    inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Uploading "${escapeHtml(file.name)}"\u2026</p>`;

    const form = new FormData();
    form.append("file", file);

    let res, data;
    try {
        res = await fetch("/upload", { method: "POST", body: form });
        data = await res.json();
    } catch (err) {
        showToast("error", `Upload failed: ${err}`);
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Upload failed: ${escapeHtml(String(err))}</p>`;
        return;
    }

    if (res.status === 413) {
        showToast("error", data.error || "File too large.");
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(data.error || "File too large.")}</p>`;
        return;
    }

    if (!data.ok) {
        showAlert(inputAlert, "error", data.error);
        uploadedRawB64 = null;
        matchedBaseEditorCount = 0;
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
}

// Spreadsheet-mode "Load example" button. Posts a small in-memory TSV
// through the same /upload pipeline as a real file pick.
const SPREADSHEET_EXAMPLE_TSV = [
    "chromosome\tposition\ttarget_base\tdesired_base",
    "12\t25245350\tT\tC",
    "12\t25245347\tT\tC",
    "17\t7675088\tT\tC",
].join("\n");

if (loadSpreadsheetExampleBtn) {
    loadSpreadsheetExampleBtn.addEventListener("click", async () => {
        if (currentMode !== "spreadsheet") setMode("spreadsheet");
        const file = new File(
            [SPREADSHEET_EXAMPLE_TSV],
            "example_targets.tsv",
            { type: "text/tab-separated-values" }
        );
        await uploadTargetsFile(file);
    });
}

// Run pipeline
runBtn.addEventListener("click", async () => {
    if (!uploadedRawB64) return;
    if (selectedEditors.size === 0) {
        showToast("error", "Select at least one base editor.");
        return;
    }
    if (currentMode === "sequence" && !customGenomeB64) {
        showToast("error", "Parse sequences before running.");
        return;
    }

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
    setSummaryPlaceholder(true, "Waiting for pipeline to complete...");
    if (resultsGuideStatus) {
        resultsGuideStatus.textContent = "";
        resultsGuideStatus.style.display = "none";
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
    form.append("base_editor", baseEditorSelectionValue());
    form.append("mode", currentMode);
    if (currentMode === "sequence") {
        form.append("custom_genome_b64", customGenomeB64 || "");
    }

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

// Render the sequence-mode preview table. Each row shows the chromosome,
// 1-based edit position, target/desired bases, and a snippet of the
// bracketed sequence so the user can spot-check the edit location.
function renderSequencePreview(rows) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);padding:16px;">No sequences to display.</p>`;
    }
    const labels = {
        chromosome: "Chromosome",
        position: "Position",
        target_base: "From",
        desired_base: "To",
        sequence_preview: "Sequence (edit highlighted)",
    };
    const cols = Object.keys(labels);
    let html = `<table class="data-table"><thead><tr>`;
    cols.forEach(c => { html += `<th>${escapeHtml(labels[c])}</th>`; });
    html += `</tr></thead><tbody>`;
    rows.forEach(row => {
        html += `<tr>`;
        cols.forEach(c => {
            let v = row[c] !== undefined && row[c] !== null ? String(row[c]) : "";
            if (c === "sequence_preview") {
                // Highlight the [X] bracket from the backend.
                const marked = escapeHtml(v).replace(
                    /\[([ACGT])\]/,
                    `<span class="sequence-edit-mark">[$1]</span>`
                );
                html += `<td><code>${marked}</code></td>`;
            } else {
                html += `<td>${escapeHtml(v)}</td>`;
            }
        });
        html += `</tr>`;
    });
    html += `</tbody></table>`;
    return html;
}

function displayColumnName(name) {
    return name.replace(/_/g, " ");
}

function renderBaseEditorTable(rows) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);">No base editor profiles available.</p>`;
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

function activeProfilesSource() {
    // With merged behavior, the Base Editors tab always shows the combined
    // list when custom profiles are loaded, otherwise just the defaults.
    if (customProfilesLoaded) {
        return "combined";
    }
    return "default";
}

async function refreshProfilesStatus(opts = {}) {
    let data;
    try {
        const res = await fetch("/base-editor-profiles/status");
        data = await res.json();
    } catch (err) {
        profilesAlert.style.display = "block";
        profilesAlert.className = "alert error";
        profilesAlert.textContent = `Failed to fetch profiles status: ${err}`;
        return;
    }

    if (!data.ok) return;

    const wasLoaded = customProfilesLoaded;
    customEditorKey = data.key || customEditorKey;
    customEditorLabel = data.label || customEditorLabel;
    customProfilesLoaded = !!data.loaded;
    customProfilesCount = (data.names || []).length;

    // When the custom-profiles set changes (upload or clear), re-fetch the
    // combined editor list so the multi-select reflects the new state.
    if (wasLoaded !== customProfilesLoaded || customProfilesLoaded) {
        try {
            await loadAvailableEditors();
        } catch (err) {
            console.warn("Failed to refresh available editors:", err);
        }
    }
}

async function refreshProfilesTable() {
    const source = activeProfilesSource();
    if (lastProfilesTableSource === source) return;
    let data;
    try {
        const res = await fetch(
            `/base-editor-table?source=${encodeURIComponent(source)}`
        );
        data = await res.json();
    } catch (err) {
        profilesTableContainer.innerHTML =
            `<p style="color:var(--danger);">Failed to load base editor profiles: ${escapeHtml(String(err))}</p>`;
        return;
    }
    if (data.ok) {
        profilesTableContainer.innerHTML = renderBaseEditorTable(data.rows);
    } else {
        profilesTableContainer.innerHTML =
            `<p style="color:var(--danger);">Failed to load base editor profiles: ${escapeHtml(data.error || "unknown error")}</p>`;
    }
    profilesSource.textContent =
        source === "combined"
            ? `Source: default + custom (${data.rows.length} profiles)`
            : `Source: default (${data.rows.length} profiles)`;
    lastProfilesTableSource = source;
}

function showProfilesAlert(kind, msg) {
    profilesAlert.style.display = "block";
    profilesAlert.className = "alert " + kind;
    profilesAlert.textContent = msg;
}

function hideProfilesAlert() {
    profilesAlert.style.display = "none";
    profilesAlert.textContent = "";
}

if (profilesUpload) {
    profilesUpload.addEventListener("change", async () => {
        const file = profilesUpload.files[0];
        if (!file) return;
        const form = new FormData();
        form.append("file", file);

        let res, data;
        try {
            res = await fetch("/base-editor-profiles/upload", {
                method: "POST",
                body: form,
            });
            data = await res.json();
        } catch (err) {
            showProfilesAlert("error", `Upload failed: ${err}`);
            profilesUpload.value = "";
            return;
        }

        if (res.status === 413) {
            showProfilesAlert("error", data.error || "File too large.");
            profilesUpload.value = "";
            return;
        }

        if (!data.ok) {
            showProfilesAlert("error", data.error);
            profilesUpload.value = "";
            return;
        }

        hideProfilesAlert();
        showProfilesAlert(
            "success",
            `Loaded ${data.count} profile(s). Switched to "${customEditorLabel}" mode.`
        );

        lastProfilesTableSource = null;
        await refreshProfilesStatus({ keepDropdownValue: false });
        await refreshProfilesTable();

        if (uploadedRawB64) {
            await runCompatibilityCheck();
        }

        profilesUpload.value = "";
        switchTab("baseEditorTab");
    });
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
        setSummaryPlaceholder(true, data.message || "No results returned.");
        if (resultsGuideStatus) {
            resultsGuideStatus.textContent = "";
            resultsGuideStatus.style.display = "none";
        }
        showToast("info", data.message || "No results returned.");
        return;
    }

    currentRecords = data.records;
    if (resultsPlaceholder) resultsPlaceholder.style.display = "none";
    resultsToolbar.style.display = "flex";
    if (resultsGuideStatus) {
        const foundTargets = new Set(
            currentRecords.map(r => `${r.chromosome ?? ""}__${r.position ?? ""}`)
        ).size;
        resultsGuideStatus.innerHTML =
            `Found guides for ` +
            `<span class="results-guide-count">${foundTargets}/${matchedBaseEditorCount}</span> ` +
            `targets with the selected base editors`;
        resultsGuideStatus.style.display = "block";
    }
    renderResultsTable();

    // Auto-switch to results tab + Guides sub-tab
    tabButtons.forEach(b => b.classList.remove("active"));
    tabPanels.forEach(p => p.classList.remove("active"));
    document.querySelector('[data-tab="resultsTab"]').classList.add("active");
    document.getElementById("resultsTab").classList.add("active");
    switchSubTab("guidesSubtab");
}

// Map the Guide score filter dropdown value to a numeric on-target score
// threshold. "all" disables the score filter (0.0 keeps every guide); the
// "okay" / "good" presets match the legacy slider defaults.
const SCORE_FILTER_THRESHOLDS = { all: 0.0, okay: 0.5, good: 0.75 };

// Dynamic filtering for the Results tab. Mirrors beguided.filtering.filters.
function getActiveFilters() {
    const key = scoreFilter.value || "all";
    const score = SCORE_FILTER_THRESHOLDS[key] ?? 0.0;
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

// Toggle the Summary sub-tab placeholder text.
function setSummaryPlaceholder(visible, text) {
    const ph = document.getElementById("summaryPlaceholder");
    if (!ph) return;
    if (visible) {
        ph.textContent = text || "Summary will appear here after the pipeline completes.";
        ph.style.display = "block";
    } else {
        ph.style.display = "none";
    }
}

// Pivoted per-base-editor summary of distinct targets visible under the current
// filters. Two rows: the top header row holds the base_editor names with
// "all" first, the second row holds the corresponding count of distinct
// (chromosome, position) targets.
function renderResultsSummary(filtered) {
    if (!resultsSummary) return;
    if (!filtered || filtered.length === 0) {
        resultsSummary.innerHTML = "";
        resultsSummary.style.display = "none";
        setSummaryPlaceholder(true, "No guides match the current filters.");
        return;
    }

    const perBaseEditor = new Map();
    const allTargets = new Set();
    for (const row of filtered) {
        const baseEditor = (row.base_editor ?? "").toString() || "(unknown)";
        const key = `${row.chromosome ?? ""}__${row.position ?? ""}`;
        if (!perBaseEditor.has(baseEditor)) perBaseEditor.set(baseEditor, new Set());
        perBaseEditor.get(baseEditor).add(key);
        allTargets.add(key);
    }

    const baseEditors = [...perBaseEditor.keys()].sort((a, b) => {
        if (a === "(unknown)") return 1;
        if (b === "(unknown)") return -1;
        return a.localeCompare(b);
    });

    let html = `<table class="data-table summary-table"><thead><tr>`;
    html += `<th>${escapeHtml(displayColumnName("base_editor"))}</th><th>all</th>`;
    for (const be of baseEditors) {
        html += `<th>${escapeHtml(be)}</th>`;
    }
    html += `</tr></thead><tbody>`;
    html += `<tr><th scope="row">targets</th>`;
    html += `<td>${allTargets.size}</td>`;
    for (const be of baseEditors) {
        html += `<td>${perBaseEditor.get(be).size}</td>`;
    }
    html += `</tr></tbody></table>`;

    resultsSummary.innerHTML = html;
    resultsSummary.style.display = "block";
    setSummaryPlaceholder(false);
}

// Filter event wiring.
bystanderFilter.addEventListener("change", () => {
    if (currentRecords.length > 0) renderResultsTable();
});
scoreFilter.addEventListener("change", () => {
    if (currentRecords.length > 0) renderResultsTable();
});
