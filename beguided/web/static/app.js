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
const gcFilter = document.getElementById("gcFilter");
const homopolymerFilter = document.getElementById("homopolymerFilter");
const trinucFilter = document.getElementById("trinucFilter");
const condenseToggle = document.getElementById("condenseToggle");
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
const loadSequenceExampleBtn = document.getElementById("loadSequenceExampleBtn");
const loadSpreadsheetExampleBtn = document.getElementById("loadSpreadsheetExampleBtn");
const profilesUpload = document.getElementById("profilesUpload");
const profilesAlert = document.getElementById("profilesAlert");
const profilesSource = document.getElementById("profilesSource");
const profilesTableContainer = document.getElementById("profilesTableContainer");

// Toast region announces status changes to screen readers.
if (toastContainer) {
    toastContainer.setAttribute("role", "status");
    toastContainer.setAttribute("aria-live", "polite");
}

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
// When true, rows that differ only in `base_editor` are collapsed into a
// single row whose `base_editor` cell lists every editor newline-separated.
// Mirrors the `condense` query param on /download/results so the downloaded
// file matches what the user sees on screen. Initialized from the checkbox
// in the results toolbar (default: checked).
let condenseEnabled = !!(condenseToggle && condenseToggle.checked);
let matchedBaseEditorCount = 0;
let customEditorKey = "custom";
let customEditorLabel = "Custom (uploaded)";
let customProfilesLoaded = false;
let lastProfilesTableSource = null;
// localStorage key for the user's last successfully uploaded custom
// profiles spreadsheet. The value is the file's bytes encoded as
// base64 so binary XLSX payloads survive intact. Bumped from
// ``beguidedCustomProfiles`` (legacy TSV text cache) to ``.v2`` when
// the upload format changed from TSV to XLSX; ``LEGACY_PROFILES_LS_KEY``
// is read once on hydrate to clean up any leftover TSV entry.
const CUSTOM_PROFILES_LS_KEY = "beguidedCustomProfiles.v2";
const LEGACY_PROFILES_LS_KEY = "beguidedCustomProfiles";
// Editor name -> editing type ("ABE" / "CBE" / ...) lookup, populated from
// /base-editor-table?source=combined. Used by the Summary bar chart to color
// per-editor bars by editing type.
const baseEditorEditingTypes = new Map();
// Live Chart.js instance for the Summary tab. Tracked so we can destroy the
// previous chart before re-rendering to avoid leaks and stacked canvases.
let summaryChartInstance = null;

// Honour the OS reduced-motion preference for Chart.js. Without this, the
// summary bar chart animates every render which can trigger vestibular
// discomfort. CSS transitions are already suppressed by the
// ``prefers-reduced-motion`` block in style.css.
(function applyReducedMotionToChart() {
    try {
        const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
        if (reduced && typeof Chart !== "undefined" && Chart && Chart.defaults) {
            Chart.defaults.animation = false;
        }
    } catch (_) {
        // matchMedia is not available (very old browsers); Chart.js falls
        // back to its default animation, which is acceptable.
    }
})();

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

// True when the editor's Cas variant is SpRY (case-insensitive). Used by the
// "All without SpRY" quick toggle and the default selection.
const isSpryCas = (cas) => /spry/i.test(cas || "");

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
        await hydrateCustomProfilesFromLocalStorage();
    } catch (err) {
        console.warn("Failed to rehydrate custom profiles:", err);
    }
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
    applyFilterTooltips();
    setMode(currentMode);
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
        origin: r.origin,
    }));
    baseEditorEditingTypes.clear();
    for (const r of data.rows) {
        if (r && r.name) baseEditorEditingTypes.set(r.name, r.editing || "");
    }
    // Default: every editor except SpRY-based ones selected.
    selectedEditors = new Set(
        availableEditors.filter(ed => !isSpryCas(ed.cas)).map(ed => ed.name)
    );
    renderCheckboxGroups();
    updateMultiselectLabel();
    updateQuickToggleVisibility();
}

// Show or hide quick-toggle buttons whose relevance depends on which
// editors are actually present: "All Custom" hides when no custom profiles
// are loaded, "All without SpRY" hides when no SpRY editors are available
// (otherwise it would be a no-op identical to "All").
function updateQuickToggleVisibility() {
    const panel = beMultiselectPanel;
    const hasCustom = availableEditors.some(ed => ed.origin === "custom");
    const customBtn = panel.querySelector(
        '.be-quick-toggles button[data-toggle="custom"]'
    );
    if (customBtn) customBtn.style.display = hasCustom ? "" : "none";

    const hasSpry = availableEditors.some(ed => isSpryCas(ed.cas));
    const spryBtn = panel.querySelector(
        '.be-quick-toggles button[data-toggle="no_spry"]'
    );
    if (spryBtn) spryBtn.style.display = hasSpry ? "" : "none";
}

// Render the checkbox groups. Order:
//   Custom (with optional Cas subgroups) -> ABE -> CBE.
function renderCheckboxGroups() {
    if (availableEditors.length === 0) {
        beCheckboxGroups.innerHTML = "";
        return;
    }

    const customEditors = availableEditors.filter(ed => ed.origin === "custom");
    const defaultEditors = availableEditors.filter(ed => ed.origin !== "custom");

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
        // Detect whether the current selection exactly matches the
        // "All without SpRY" preset (every non-SpRY editor selected, no
        // SpRY editor selected) and label accordingly. Other partial
        // selections fall through to the generic "N selected" label.
        const nonSpryEditors = availableEditors.filter(ed => !isSpryCas(ed.cas));
        const nonSpryTotal = nonSpryEditors.length;
        const nonSpryNames = new Set(nonSpryEditors.map(ed => ed.name));
        const selectedIsNonSpry = nonSpryTotal > 0
            && [...selectedEditors].every(name => nonSpryNames.has(name))
            && n === nonSpryTotal;
        if (selectedIsNonSpry) {
            beMultiselectLabel.textContent = `All without SpRY (${nonSpryTotal})`;
            return;
        }
        beMultiselectLabel.textContent = `${n} selected`;
    }
}

// Quick toggle: all / no_spry / abe / cbe / custom / none.
function applyQuickToggle(toggle) {
    if (toggle === "all") {
        selectedEditors = new Set(availableEditors.map(r => r.name));
    } else if (toggle === "no_spry") {
        selectedEditors = new Set(
            availableEditors.filter(r => !isSpryCas(r.cas)).map(r => r.name)
        );
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
                .filter(r => r.origin === "custom")
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
        if (ev.key === "Escape" && !beMultiselectPanel.hasAttribute("hidden")) {
            closeMultiselectPanel();
            beMultiselectBtn.focus();
            return;
        }
        if (beMultiselectPanel.hasAttribute("hidden")) return;
        if (ev.key === "Home" || ev.key === "End") {
            const boxes = [...beCheckboxGroups.querySelectorAll('input[type="checkbox"]')];
            if (boxes.length === 0) return;
            ev.preventDefault();
            const target = ev.key === "Home" ? boxes[0] : boxes[boxes.length - 1];
            target.focus();
        } else if (ev.key === "Tab") {
            // Focus trap: keep Tab cycling within the open panel.
            const focusables = [
                ...beMultiselectPanel.querySelectorAll(
                    'button, input[type="checkbox"]'
                ),
            ].filter(el => !el.disabled && el.offsetParent !== null);
            if (focusables.length === 0) return;
            const first = focusables[0];
            const last = focusables[focusables.length - 1];
            const active = document.activeElement;
            if (ev.shiftKey && active === first) {
                ev.preventDefault();
                last.focus();
            } else if (!ev.shiftKey && active === last) {
                ev.preventDefault();
                first.focus();
            }
        }
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
        const isActive = b.dataset.tab === tabId;
        b.classList.toggle("active", isActive);
        b.setAttribute("aria-selected", isActive ? "true" : "false");
        b.tabIndex = isActive ? 0 : -1;
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
        const isActive = b.dataset.subtab === subtabId;
        b.classList.toggle("active", isActive);
        b.setAttribute("aria-selected", isActive ? "true" : "false");
        b.tabIndex = isActive ? 0 : -1;
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
    modeSpreadsheetBtn.setAttribute("aria-selected", mode === "spreadsheet" ? "true" : "false");
    modeSequenceBtn.setAttribute("aria-selected", mode === "sequence" ? "true" : "false");
    modeSpreadsheetBtn.tabIndex = mode === "spreadsheet" ? 0 : -1;
    modeSequenceBtn.tabIndex = mode === "sequence" ? 0 : -1;
    spreadsheetPanel.style.display = mode === "spreadsheet" ? "" : "none";
    sequencePanel.style.display = mode === "sequence" ? "" : "none";

    // Move the shared Run button into the active input panel's action row so
    // it stays visually anchored next to the input controls. The 'right' class
    // pushes it to the far edge of the row in spreadsheet mode; in sequence
    // mode it sits next to the Load example button, so we strip it.
    const targetRow = mode === "spreadsheet"
        ? spreadsheetPanel.querySelector(".input-actions-row")
        : sequencePanel.querySelector(".input-meta-row");
    if (targetRow) {
        targetRow.appendChild(runBtn);
    }
    runBtn.classList.toggle("right", mode === "spreadsheet");

    // Reset state when switching modes: a previously uploaded file is
    // irrelevant for sequence input and vice-versa. In sequence mode the Run
    // button starts enabled because clicking it triggers the parse + check +
    // run flow itself; in spreadsheet mode it stays disabled until a file is
    // uploaded (handled by uploadTargetsFile → runCompatibilityCheck).
    uploadedRawB64 = null;
    customGenomeB64 = null;
    matchedBaseEditorCount = 0;
    runBtn.disabled = (mode === "spreadsheet");
    compatibilityInfo.textContent = "";
    checkAlert.style.display = "none";
    inputAlert.style.display = "none";
    inputPanel.innerHTML = `<p style="color: var(--muted);padding:16px;">${
        mode === "spreadsheet"
            ? "Upload a targets file to preview input data."
            : "Paste sequences above and click Run to parse and start the pipeline."
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

// Sequence-mode parse helper. Posts the textarea content to /sequence/parse
// and returns a structured result. No UI side effects — callers decide how
// to surface errors and whether to switch tabs.
async function parseSequenceText() {
    const text = sequenceInput.value;
    if (!text.trim()) {
        return { ok: false, error: "Paste at least one sequence." };
    }
    const form = new FormData();
    form.append("text", text);
    let res, data;
    try {
        res = await fetch("/sequence/parse", { method: "POST", body: form });
        data = await res.json();
    } catch (err) {
        return { ok: false, error: `Parse failed: ${err}` };
    }
    if (res.status === 413 || !data.ok) {
        return { ok: false, error: data.error || "Parse failed." };
    }
    return {
        ok: true,
        raw_b64: data.raw_b64,
        custom_genome_b64: data.custom_genome_b64,
        preview: data.preview,
        row_count: data.row_count,
    };
}

// Sequence-mode "Load example" button. Populates the textarea with a small
// multi-record FASTA (matches the example shown in the textarea's explanation
// text). Parsing now happens on Run click, not here.
const SEQUENCE_EXAMPLE_FASTA =
    ">my_target_1\n" +
    "ATAGAGTATAGACCC[A>G]CTATCCAGTCTACCTAGGGAATC\n" +
    ">my_target_2\n" +
    "GGGCATAGCCATGTAGCATT[T>C]GCAATTGGCCATGCATTTACC";

if (loadSequenceExampleBtn) {
    loadSequenceExampleBtn.addEventListener("click", () => {
        if (currentMode !== "sequence") setMode("sequence");
        sequenceInput.value = SEQUENCE_EXAMPLE_FASTA;
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
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(data.error || "Upload failed.")}</p>`;
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
    if (selectedEditors.size === 0) {
        showToast("error", "Select at least one base editor.");
        return;
    }

    // Sequence mode: parse on Run if not already parsed. On parse error, route
    // the user to the Input tab so the error is visible; on success, render
    // the preview silently into inputPanel and continue without switching tabs.
    let parsedThisClick = false;
    if (currentMode === "sequence" && !customGenomeB64) {
        inputAlert.style.display = "none";
        checkAlert.style.display = "none";
        compatibilityInfo.textContent = "";
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Parsing sequences\u2026</p>`;
        const parsed = await parseSequenceText();
        if (!parsed.ok) {
            switchTab("inputTab");
            showAlert(inputAlert, "error", parsed.error);
            inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(parsed.error)}</p>`;
            runBtn.disabled = false;
            return;
        }
        uploadedRawB64 = parsed.raw_b64;
        customGenomeB64 = parsed.custom_genome_b64;
        inputPanel.innerHTML = renderSequencePreview(parsed.preview);
        parsedThisClick = true;
    }

    if (!uploadedRawB64) return;

    // When parse just happened, re-run the compatibility check before
    // submitting so the user is alerted to any target/editor mismatches. On
    // mismatch, route to the Input tab where the checkAlert is shown.
    if (parsedThisClick) {
        await runCompatibilityCheck();
        if (runBtn.disabled) {
            switchTab("inputTab");
            return;
        }
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

function renderTable(rows, columns, displayColumns = null, groupKeyFn = null, cellFormatter = null, headerTooltips = null, rowClassFn = null) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);padding:16px;">No data to display.</p>`;
    }
    const headers = displayColumns || columns;
    let html = `<table class="data-table"><thead><tr>`;
    headers.forEach((c, i) => {
        const dataCol = escapeHtml(columns[i] ?? c);
        const tip = headerTooltips ? headerTooltips[columns[i]] : null;
        const tipAttr = tip ? ` title="${escapeHtml(tip)}"` : "";
        const tipIcon = tip ? ` <span class="th-info" aria-hidden="true">\u24D8</span>` : "";
        html += `<th scope="col" data-col="${dataCol}"${tipAttr}>${escapeHtml(c)}${tipIcon}</th>`;
    });
    html += `</tr></thead><tbody>`;
    let prevKey = null;
    let stripeToggle = false;
    rows.forEach(row => {
        let cls = rowClassFn ? (rowClassFn(row) || "") : "";
        if (groupKeyFn) {
            const key = groupKeyFn(row);
            if (key !== prevKey) {
                stripeToggle = !stripeToggle;
                prevKey = key;
            }
            cls = (cls ? cls + " " : "") + (stripeToggle ? "target-a" : "target-b");
        }
        html += cls ? `<tr class="${cls}">` : `<tr>`;
        columns.forEach(c => {
            const v = row[c] !== undefined ? String(row[c]) : "";
            const formatted = cellFormatter ? cellFormatter(row, c, v) : null;
            const content = formatted !== null ? formatted : escapeHtml(v);
            if (c === "base_editor") {
                html += `<td data-col="base_editor"><div class="be-cell">${content}</div></td>`;
            } else {
                html += `<td data-col="${escapeHtml(c)}">${content}</td>`;
            }
        });
        html += `</tr>`;
    });
    html += `</tbody></table>`;
    return html;
}

// Parse the comma-separated ``bystander_positions`` string into a Set of
// 0-based sequence indices. Bystander positions are 1-based from the 5'
// end of the guide (same convention as ``editing_position``). Tokens that
// are non-numeric, out of range, or duplicates are silently dropped.
function parseBystanderIndices(raw, seqLen) {
    const idxs = new Set();
    if (!raw) return idxs;
    for (const tok of String(raw).split(",")) {
        const n = Number(tok.trim());
        if (Number.isInteger(n) && n >= 1 && n <= seqLen) {
            idxs.add(n - 1);
        }
    }
    return idxs;
}

// Format DNA-sequence cells in the results table so they share a uniform
// monospace look. The guide RNA sequence gets the base at the 1-based
// ``editing_position`` (from the 5' end of the guide) highlighted in blue
// and any bases listed in ``bystander_positions`` highlighted in dark red;
// the PAM sequence gets the same monospace treatment but no highlight
// since it has no positional metadata.
function formatDnaSequenceCell(row, column, value) {
    if (column !== "guide_rna_sequence" && column !== "pam_sequence") {
        return null;
    }
    const seq = String(value ?? "");
    const safeSeq = escapeHtml(seq);
    if (column === "pam_sequence") {
        return `<code>${safeSeq}</code>`;
    }
    const editPos = Number(row.editing_position);
    const editIdx = Number.isInteger(editPos) ? editPos - 1 : -1;
    const bystanderIdxs = parseBystanderIndices(row.bystander_positions, seq.length);
    if (editIdx < 0 && bystanderIdxs.size === 0) {
        return `<code>${safeSeq}</code>`;
    }
    const parts = [];
    for (let i = 0; i < seq.length; i++) {
        if (i === editIdx) {
            parts.push(`<span class="guide-edit-base" aria-label="edited base ${safeSeq[i]} at position ${i + 1}">${safeSeq[i]}</span>`);
        } else if (bystanderIdxs.has(i)) {
            parts.push(`<span class="guide-bystander-base" aria-label="bystander base ${safeSeq[i]} at position ${i + 1}">${safeSeq[i]}</span>`);
        } else {
            parts.push(safeSeq[i]);
        }
    }
    return `<code>${parts.join("")}</code>`;
}

function escapeHtml(str) {
    const div = document.createElement("div");
    div.textContent = str;
    return div.innerHTML;
}

// Render the guidefinder-derived columns. ``gc_content`` is shown as a
// percentage with one decimal (the column value is percent 0-100); a missing
// value renders as an em-dash. ``homopolymer`` is rendered as a small
// coloured pill ("poly" or "ok") so flagged guides stand out at a glance.
// ``trinuc_edit_score`` (Marquart et al. 2021, ABE8e, 0..1) is shown as a
// 2-decimal number; ``NaN``/null (CBE rows) render as an em-dash.
function formatSequencePropsCell(row, column, value) {
    if (column === "gc_content") {
        if (value === undefined || value === null || value === "") {
            return `<span class="muted">&mdash;</span>`;
        }
        const num = Number(value);
        if (!Number.isFinite(num)) {
            return escapeHtml(String(value));
        }
        return `${num.toFixed(1)}%`;
    }
    if (column === "homopolymer") {
        if (value === true || value === "true") {
            return `<span class="badge badge-warn">poly</span>`;
        }
        if (value === false || value === "false") {
            return `<span class="muted">ok</span>`;
        }
        return `<span class="muted">&mdash;</span>`;
    }
    if (column === "trinuc_edit_score") {
        if (value === undefined || value === null || value === "") {
            return `<span class="muted">&mdash;</span>`;
        }
        const num = Number(value);
        if (!Number.isFinite(num)) {
            // CBE rows have NaN; missing values render the same way.
            return `<span class="muted">&mdash;</span>`;
        }
        return num.toFixed(2);
    }
    return null;
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
    cols.forEach(c => { html += `<th scope="col">${escapeHtml(labels[c])}</th>`; });
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

const SEQUENCE_MODE_COLUMN_LABELS = {
    chromosome: "sequence",
};

// Display-label overrides for result columns, applied in all modes.
// The underlying data column name stays unchanged (e.g. ``trinuc_edit_score``);
// only the header text shown to the user is remapped.
const COLUMN_LABELS = {
    trinuc_edit_score: "context edit score",
};

// Tooltip text shown on hover for results-table headers. Keyed by the
// snake_case column name (as emitted by the pipeline), so the lookup works
// regardless of the display label. Only the results-table call site passes
// this map to renderTable(); the upload-preview table is left untouched.
const RESULTS_COLUMN_TOOLTIPS = {
    gene_edit_desired:
        "Desired edit at the protein level, formatted as '<gene> <Target_AA><codon_number><Desired_AA>' (e.g. 'TP53 *306R'). Codon numbering follows the MANE Select CDS frame; empty when the target is not in a CDS or no gene name resolves.",
    editing_position:
        "1-based position of the edited base, counted from the 5' end of the guide RNA (not the genomic coordinate).",
    bystander_positions:
        "Comma-separated 1-based positions from the 5' end of the guide where additional bases may be edited as bystanders.",
    bystander_effect:
        "Predicted coding-sequence consequence of the bystander edits (e.g. silent / missense / nonsense), assessed against the MANE Select transcript codon frame. Empty when no GTF was provided.",
    guide_score:
        "CRISPRater on-target activity score (0-1). Derived from spCas9 training data and used here as a general guide-quality proxy.",
    gc_content:
        "GC fraction of the guide RNA spacer, expressed as a percentage. Guides with very low or very high GC tend to bind less efficiently.",
    homopolymer:
        "True when the spacer contains a run of >=4 identical bases (e.g. AAAA). Homopolymer runs can hurt guide specificity and efficiency.",
    trinuc_edit_score:
        "Marquart et al. 2021 normalised ABE8e editing proportion for the local X-A-Y trinucleotide context (the bases immediately before and after the target A). 0-1; NaN for CBE guides.",
};

// Maps each results-toolbar filter <select>/<input> id to the column whose
// RESULTS_COLUMN_TOOLTIPS entry describes it. Tooltips are wired up at init
// so the table header and the filter label always share the same text.
const FILTER_TOOLTIP_COLUMN = {
    bystanderFilter: "bystander_effect",
    scoreFilter: "guide_score",
    gcFilter: "gc_content",
    homopolymerFilter: "homopolymer",
    trinucFilter: "trinuc_edit_score",
};

// Apply RESULTS_COLUMN_TOOLTIPS entries to the corresponding filter <label>s
// in the results toolbar. Sets both the native `title` attribute (browser
// tooltip on hover) and appends the same info glyph used by the table headers
// so users can see at a glance which filters have help text.
function applyFilterTooltips() {
    for (const [id, col] of Object.entries(FILTER_TOOLTIP_COLUMN)) {
        const tip = RESULTS_COLUMN_TOOLTIPS[col];
        if (!tip) continue;
        const label = document.querySelector(`label[for="${id}"]`);
        if (!label) continue;
        label.title = tip;
        if (!label.querySelector(".th-info")) {
            const glyph = document.createElement("span");
            glyph.className = "th-info";
            glyph.setAttribute("aria-hidden", "true");
            glyph.textContent = "\u24D8";
            label.appendChild(glyph);
        }
    }
}

function columnLabel(name) {
    if (currentMode === "sequence"
        && Object.prototype.hasOwnProperty.call(SEQUENCE_MODE_COLUMN_LABELS, name)) {
        return SEQUENCE_MODE_COLUMN_LABELS[name];
    }
    if (Object.prototype.hasOwnProperty.call(COLUMN_LABELS, name)) {
        return COLUMN_LABELS[name];
    }
    return displayColumnName(name);
}

// Reorder columns to put the specified column first, preserving order of others
function reorderColumns(columns, firstColumn) {
    const index = columns.indexOf(firstColumn);
    if (index === -1) return columns;
    const reordered = [...columns];
    reordered.splice(index, 1);  // Remove from current position
    reordered.unshift(firstColumn);  // Add to front
    return reordered;
}

function displayColumnName(name) {
    return name.replace(/_/g, " ");
}

// Render the ``source`` cell as a clickable link when the backend supplies
// a parsed URL; otherwise fall back to the raw (HTML-escaped) text. The
// backend only ever emits http/https URLs here, but we still escape both
// the href and the link text so a stray value cannot inject markup.
function formatSourceCell(value) {
    const src = (value && typeof value === "object") ? value : { text: "", url: null };
    const text = src.text == null ? "" : String(src.text);
    if (src.url) {
        return `<a href="${escapeHtml(src.url)}" target="_blank" rel="noopener noreferrer">${escapeHtml(text)}</a>`;
    }
    return escapeHtml(text);
}

function renderBaseEditorTable(rows) {
    if (!rows || rows.length === 0) {
        return `<p style="color:var(--muted);">No base editor profiles available.</p>`;
    }
    const cols = ["name", "editing", "general", "high", "pam", "cas", "source"];
    const labels = {
        name: "Name",
        editing: "Editing",
        general: "General window",
        high: "High-activity window",
        pam: "PAM",
        cas: "Cas",
        source: "Source",
    };
    let html = `<table class="data-table"><thead><tr>`;
    cols.forEach(c => {
        html += `<th scope="col">${escapeHtml(labels[c])}</th>`;
    });
    html += `</tr></thead><tbody>`;
    rows.forEach(row => {
        html += `<tr>`;
        cols.forEach(c => {
            if (c === "source") {
                html += `<td>${formatSourceCell(row[c])}</td>`;
            } else {
                const v = row[c] !== undefined ? String(row[c]) : "";
                html += `<td>${escapeHtml(v)}</td>`;
            }
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

function readCustomProfilesFromLocalStorage() {
    try {
        return localStorage.getItem(CUSTOM_PROFILES_LS_KEY);
    } catch (_) {
        return null;
    }
}

// Encode an ArrayBuffer of arbitrary profile bytes (XLSX is binary) as
// a base64 string so localStorage can carry it. UTF-8 round-trip via
// String.fromCharCode / atob would corrupt non-ASCII bytes, which a
// base64 detour avoids.
function arrayBufferToBase64(buffer) {
    const bytes = new Uint8Array(buffer);
    let binary = "";
    const chunk = 0x8000;
    for (let i = 0; i < bytes.length; i += chunk) {
        binary += String.fromCharCode.apply(
            null, bytes.subarray(i, i + chunk)
        );
    }
    return btoa(binary);
}

function base64ToArrayBuffer(b64) {
    const binary = atob(b64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) {
        bytes[i] = binary.charCodeAt(i);
    }
    return bytes.buffer;
}

function writeCustomProfilesToLocalStorage(b64) {
    try {
        localStorage.setItem(CUSTOM_PROFILES_LS_KEY, b64);
        return true;
    } catch (_) {
        return false;
    }
}

function clearCustomProfilesFromLocalStorage() {
    try {
        localStorage.removeItem(CUSTOM_PROFILES_LS_KEY);
    } catch (_) {
    }
}

async function hydrateCustomProfilesFromLocalStorage() {
    // One-shot migration: drop the legacy TSV-text cache so we don't
    // keep re-uploading a payload the server no longer accepts.
    try {
        localStorage.removeItem(LEGACY_PROFILES_LS_KEY);
    } catch (_) {
    }
    const b64 = readCustomProfilesFromLocalStorage();
    if (!b64) return;
    try {
        const buf = base64ToArrayBuffer(b64);
        const form = new FormData();
        form.append(
            "file",
            new Blob([buf], {
                type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            }),
            "custom_profiles.xlsx"
        );
        const res = await fetch("/base-editor-profiles/upload", {
            method: "POST",
            body: form,
        });
        const data = await res.json().catch(() => ({ ok: false }));
        if (!data || !data.ok) {
            clearCustomProfilesFromLocalStorage();
        }
    } catch (err) {
        console.warn("Failed to rehydrate custom profiles from localStorage:", err);
    }
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
        const fileBuffer = await file.arrayBuffer();
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

        writeCustomProfilesToLocalStorage(arrayBufferToBase64(fileBuffer));

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
                if (document.activeElement === cancelBtn) {
                    runBtn.focus();
                }
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

function targetKey(r) {
    const chrom = r.chromosome ?? r.sequence ?? "";
    const pos = r.position ?? "";
    return `${chrom}__${pos}`;
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
            currentRecords.filter(r => !isNoGuideRow(r)).map(targetKey)
        ).size;
        resultsGuideStatus.innerHTML =
            foundTargets > 0
                ? `Found guides for ` +
                  `<span class="results-guide-count">${foundTargets}/${matchedBaseEditorCount}</span> ` +
                  `targets with the selected base editors`
                : `No guides found for the selected base editors`;
        resultsGuideStatus.style.display = "block";
    }
    renderResultsTable();

    // Auto-switch to results tab + Guides sub-tab
    switchTab("resultsTab");
    switchSubTab("guidesSubtab");
    if (resultsGuideStatus && resultsGuideStatus.style.display !== "none") {
        resultsGuideStatus.setAttribute("tabindex", "-1");
        resultsGuideStatus.focus();
    }
}

// Map the Guide score filter dropdown value to a numeric on-target score
// threshold. "all" disables the score filter (0.0 keeps every guide); the
// "acceptable" / "good" presets match the legacy slider defaults.
const SCORE_FILTER_THRESHOLDS = { all: 0.0, acceptable: 0.5, good: 0.75 };

// Map the Context edit score filter dropdown value to a numeric threshold.
// Mirrors the trinucleotide edit-context score (Marquart et al. 2021, ABE8e,
// 0..1). "all" disables the filter; CBE rows have a null/NaN score and are
// always kept, mirroring the server's ``df.isna() | (df >= threshold)`` mask.
const CONTEXT_EDIT_SCORE_THRESHOLDS = {
    all: 0.0,
    acceptable: 0.25,
    good: 0.5,
    very_good: 0.75,
};

// GC content tier ranges (percent 0-100). "all" disables the filter;
// "acceptable" keeps guides with GC in 30-70%, "good" narrows to 40-60%.
// The frontend uses "all" for consistency with the other dropdown filters;
// it's translated to "off" before being sent to the server.
const GC_FILTER_RANGES = {
    all: null,
    good: [40, 60],
    acceptable: [30, 70],
};

// Dynamic filtering for the Results tab. Mirrors beguided.filtering.filters.
function getActiveFilters() {
    const key = scoreFilter.value || "all";
    const score = SCORE_FILTER_THRESHOLDS[key] ?? 0.0;
    // Context edit-score filter (Marquart et al. 2021, ABE8e, 0..1).
    // "all" disables the filter (0.0 keeps every guide). CBE rows have a
    // null/NaN ``trinuc_edit_score`` and are always kept, mirroring the
    // server's ``df.isna() | (df >= threshold)`` mask.
    const trinucKey = (trinucFilter && trinucFilter.value) || "all";
    const trinuc = CONTEXT_EDIT_SCORE_THRESHOLDS[trinucKey] ?? 0.0;
    return {
        mode: bystanderFilter.value,
        score,
        gc: (gcFilter && gcFilter.value) || "all",
        homopolymer: !!(homopolymerFilter && homopolymerFilter.checked),
        trinuc,
    };
}

// Collapse records that differ only in `base_editor` into one row whose
// `base_editor` cell lists every editor newline-separated (alphabetical).
// The grouping key is every column except `base_editor`, matching the
// server-side `condense_df` helper in beguided/web/condense.py so the
// downloaded file matches what the user sees on screen.
function condenseRecords(records) {
    if (!records || records.length === 0) return records || [];
    if (!records.every(r => r && Object.prototype.hasOwnProperty.call(r, "base_editor"))) {
        return records;
    }
    const buckets = new Map();
    const order = [];
    for (const r of records) {
        const keyParts = [];
        for (const k of Object.keys(r)) {
            if (k === "base_editor") continue;
            keyParts.push(_condenseKeyValue(k, r[k]));
        }
        const key = keyParts.join("\u0001");
        if (!buckets.has(key)) {
            buckets.set(key, []);
            order.push(key);
        }
        buckets.get(key).push(r);
    }
    const out = [];
    for (const key of order) {
        const members = buckets.get(key);
        const merged = { ...members[0] };
        const editors = [...new Set(members.map(m => String(m.base_editor)))].sort();
        merged.base_editor = editors.join("\n");
        out.push(merged);
    }
    return out;
}

function _condenseKeyValue(_key, value) {
    if (value === undefined) return "__u__";
    if (value === null) return "__n__";
    if (typeof value === "object") return JSON.stringify(value);
    return String(value);
}

function filterRecords(records, mode, minScore, gcMode = "all", dropHomopolymers = false, trinucMin = null) {
    const gcRange = GC_FILTER_RANGES[gcMode] ?? null;
    return records.filter(row => {
        // On-target score filter. Missing scores are kept.
        if (minScore > 0) {
            const raw = row.guide_score;
            if (raw !== undefined && raw !== null && raw !== "") {
                const num = Number(raw);
                if (Number.isFinite(num) && num < minScore) return false;
            }
        }

        // GC content tier filter. Missing GC is kept, matching the server.
        if (gcRange) {
            const rawGc = row.gc_content;
            if (rawGc !== undefined && rawGc !== null && rawGc !== "") {
                const num = Number(rawGc);
                if (Number.isFinite(num) && (num < gcRange[0] || num > gcRange[1])) {
                    return false;
                }
            }
        }

        // Homopolymer filter. Missing flag is kept.
        if (dropHomopolymers) {
            const rawHp = row.homopolymer;
            if (rawHp === true || rawHp === "true") return false;
        }

        // Context edit-score filter (ABE only; CBE rows kept).
        // Missing / NaN values are kept, matching the server.
        if (trinucMin !== null && trinucMin > 0) {
            const rawTrinuc = row.trinuc_edit_score;
            if (rawTrinuc !== undefined && rawTrinuc !== null && rawTrinuc !== "") {
                const num = Number(rawTrinuc);
                if (Number.isFinite(num) && num < trinucMin) return false;
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

// Placeholder-row markers for targets that ended up without a guide. Mirror
// beguided.models.PLACEHOLDER_MARKERS so the client can recognise pipeline
// placeholders (added in the pipeline's very last step) and build its own
// "**Guides removed by filter**" rows when the dynamic filters below remove
// every guide of a target.
const NO_GUIDE_FOUND = "**No guide possible**";
const BASE_CHANGE_INCOMPATIBLE = "**Base change incompatible**";
const GUIDES_REMOVED_BY_FILTER = "**Guides removed by filter**";
const PLACEHOLDER_MARKERS = new Set([
    NO_GUIDE_FOUND,
    BASE_CHANGE_INCOMPATIBLE,
    GUIDES_REMOVED_BY_FILTER,
]);
// Identity/base columns copied verbatim onto a placeholder row; every other
// result column is emptied. Covers spreadsheet mode (chromosome/position,
// target_base/desired_base) and sequence mode (sequence, desired_edit, no
// position). Mirrors filters._PLACEHOLDER_IDENTITY_COLUMNS.
const PLACEHOLDER_IDENTITY_COLUMNS = new Set([
    "chromosome", "sequence", "position", "target_base", "desired_base", "desired_edit",
]);

function isNoGuideRow(row) {
    const seq = row && row.guide_rna_sequence;
    return seq !== undefined && seq !== null && PLACEHOLDER_MARKERS.has(String(seq));
}

// Build a placeholder row for a target whose guides were all removed by the
// dynamic filters, sharing the server-side semantics: identity/base columns
// come from the template row, result columns are emptied, and the marker is
// preserved when the template row was already a placeholder.
function makeNoGuideRow(template) {
    const row = {};
    for (const k of Object.keys(template)) {
        row[k] = PLACEHOLDER_IDENTITY_COLUMNS.has(k) ? template[k] : null;
    }
    const sourceSeq = template.guide_rna_sequence;
    row.guide_rna_sequence = PLACEHOLDER_MARKERS.has(String(sourceSeq))
        ? sourceSeq
        : GUIDES_REMOVED_BY_FILTER;
    return row;
}

// After filtering, append one placeholder row per target (unique
// ``targetKey``) that has no surviving row. Mirrors
// ``beguided.filtering.filters._fill_missing_target_rows`` so the on-screen
// table matches the downloaded file.
function fillMissingTargets(filtered, preFilter) {
    if (!preFilter || preFilter.length === 0) return filtered;
    const present = new Set(filtered.map(targetKey));
    // First row per target in the pre-filter set is the template for the
    // placeholder (it carries the identity/base info even when the target's
    // guides were all removed).
    const firstIndex = new Map();
    for (let i = 0; i < preFilter.length; i++) {
        const key = targetKey(preFilter[i]);
        if (!firstIndex.has(key)) firstIndex.set(key, i);
    }
    const out = filtered.slice();
    for (const [key, idx] of firstIndex) {
        if (present.has(key)) continue;
        out.push(makeNoGuideRow(preFilter[idx]));
    }
    return out;
}

function renderResultsTable() {
    const wrap = resultsPanel.querySelector(".results-table-wrap");
    if (!wrap) {
        const div = document.createElement("div");
        div.className = "results-table-wrap";
        resultsPanel.appendChild(div);
    }
    const target = resultsPanel.querySelector(".results-table-wrap");

    const { mode, score, gc, homopolymer, trinuc } = getActiveFilters();
    // Condense BEFORE filtering so the filter sees the same row shape the
    // download endpoint will return. When the toggle is off this is a no-op
    // (condenseRecords passes records through unchanged).
    const condensed = condenseEnabled ? condenseRecords(currentRecords) : currentRecords;
    const baseFiltered = filterRecords(condensed, mode, score, gc, homopolymer, trinuc);
    // Re-add one placeholder row per target whose guides were all removed by
    // the active filters (mirrors the server-side download fill).
    const filtered = fillMissingTargets(baseFiltered, currentRecords);
    // Summary is intentionally NOT affected by the condense toggle: compute a
    // separate filtered view from the uncondensed records so per-editor counts
    // always reflect one editor per row. Filters still apply.
    const summaryFiltered = filterRecords(currentRecords, mode, score, gc, homopolymer, trinuc);

    target.classList.toggle("condensed", condenseEnabled && condensed.length < currentRecords.length);

    if (filtered.length === 0) {
        target.innerHTML = `<p style="color:var(--muted);padding:16px;">No guides match the current filters.</p>`;
    } else {
        const columns = Object.keys(filtered[0]);
        // Move 'sequence' column to the front in sequence mode
        const finalColumns = currentMode === "sequence"
            ? reorderColumns(columns, "sequence")
            : columns;
        const displayColumns = finalColumns.map(columnLabel);
        target.innerHTML = renderTable(
            filtered,
            finalColumns,
            displayColumns,
            targetKey,
            (row, column, value) => {
                const seqProps = formatSequencePropsCell(row, column, value);
                if (seqProps !== null) return seqProps;
                return formatDnaSequenceCell(row, column, value);
            },
            RESULTS_COLUMN_TOOLTIPS,
            (row) => isNoGuideRow(row) ? "no-guide" : null
        );
    }

    renderResultsSummary(summaryFiltered);

    filterSummary.textContent = formatFilterSummary(
        filtered.length, condensed.length, currentRecords.length
    );

    if (currentJobId && currentRecords.length > 0) {
        const params = new URLSearchParams({
            job_id: currentJobId,
            filter_mode: mode,
            guide_score_threshold: score.toString(),
            // Translate the frontend's "all" to the server's "off" so the
            // download endpoint's VALID_GC_FILTERS check accepts the value.
            gc_content_filter: gc === "all" ? "off" : gc,
        });
        if (homopolymer) params.append("filter_homopolymers", "1");
        if (trinuc > 0) params.append("trinuc_edit_score_min", trinuc.toString());
        if (condenseEnabled) params.append("condense", "1");
        const baseHref = "/download/results?" + params.toString();
        downloadResults.href = baseHref;
        downloadResults.style.display = "inline-flex";
        downloadResultsXlsx.href = baseHref + "&format=xlsx";
        downloadResultsXlsx.style.display = "inline-flex";
    }
}

function formatFilterSummary(shown, afterCondense, total) {
    if (condenseEnabled && afterCondense < total) {
        return `Showing ${shown} of ${total} (condensed: ${total} \u2192 ${afterCondense} guides)`;
    }
    return `Showing ${shown} of ${total}`;
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

// Read a CSS custom property from :root, with a hard fallback for environments
// where the variable is missing.
function cssVar(name, fallback) {
    const v = getComputedStyle(document.documentElement).getPropertyValue(name);
    return (v && v.trim()) || fallback;
}

// Resolve the chart bar color for a given label. "all" gets the accent color,
// per-editor bars get the color of their editing type (ABE/CBE), and anything
// unrecognized falls back to the muted color.
function summaryBarColor(label) {
    if (label === "all") return cssVar("--chart-all", "#475569");
    const t = (baseEditorEditingTypes.get(label) || "").toUpperCase();
    if (t === "ABE") return cssVar("--chart-abe", "#2563eb");
    if (t === "CBE") return cssVar("--chart-cbe", "#f59e0b");
    return cssVar("--chart-unknown", "#94a3b8");
}

// Inline Chart.js plugin: render the numeric value at the right end of every
// bar so users can read exact counts without hovering. Avoids pulling in the
// separate chartjs-plugin-datalabels dependency.
const summaryValueLabelsPlugin = {
    id: "summaryValueLabels",
    afterDatasetsDraw(chart) {
        const { ctx } = chart;
        const meta = chart.getDatasetMeta(0);
        if (!meta || !meta.data) return;
        ctx.save();
        ctx.font = "600 12px system-ui, -apple-system, Segoe UI, Roboto, sans-serif";
        ctx.fillStyle = cssVar("--text", "#1e293b");
        ctx.textBaseline = "middle";
        ctx.textAlign = "left";
        const values = chart.data.datasets[0].data || [];
        meta.data.forEach((bar, i) => {
            const v = values[i];
            if (v === undefined || v === null) return;
            const label = `${v}`;
            const x = bar.x + 6;
            const y = bar.y;
            ctx.fillText(label, x, y);
        });
        ctx.restore();
    },
};

// Pivoted per-base-editor summary of distinct targets visible under the current
// filters. Rendered as a horizontal bar chart: "all" sits on top, per-editor
// bars follow in alphabetical order, with the editing type driving bar color
// and the count rendered at the right end of each bar.
function renderResultsSummary(filtered) {
    if (!resultsSummary) return;
    if (!filtered || filtered.length === 0) {
        destroySummaryChart();
        resultsSummary.innerHTML = "";
        resultsSummary.style.display = "none";
        setSummaryPlaceholder(true, "No guides match the current filters.");
        return;
    }

    // Targets without a guide RNA are intentionally excluded from the
    // per-editor summary; they are still visible in the Guides table.
    const realRows = filtered.filter(r => !isNoGuideRow(r));
    if (realRows.length === 0) {
        destroySummaryChart();
        resultsSummary.innerHTML = "";
        resultsSummary.style.display = "none";
        setSummaryPlaceholder(true, "No guides match the current filters.");
        return;
    }

    const perBaseEditor = new Map();
    const allTargets = new Set();
    for (const row of realRows) {
        const baseEditor = (row.base_editor ?? "").toString() || "(unknown)";
        const key = targetKey(row);
        if (!perBaseEditor.has(baseEditor)) perBaseEditor.set(baseEditor, new Set());
        perBaseEditor.get(baseEditor).add(key);
        allTargets.add(key);
    }

    const baseEditors = [...perBaseEditor.keys()].sort((a, b) => {
        if (a === "(unknown)") return 1;
        if (b === "(unknown)") return -1;
        return a.localeCompare(b);
    });

    const labels = ["all", ...baseEditors];
    const data = [allTargets.size, ...baseEditors.map(be => perBaseEditor.get(be).size)];
    const colors = labels.map(summaryBarColor);

    resultsSummary.innerHTML =
        `<div class="summary-chart-wrap"><canvas id="summaryChart" aria-label="Targets per base editor"></canvas></div>`;
    resultsSummary.style.display = "block";
    setSummaryPlaceholder(false);

    drawSummaryChart(labels, data, colors);
}

function destroySummaryChart() {
    if (summaryChartInstance) {
        try { summaryChartInstance.destroy(); } catch (_) { /* noop */ }
        summaryChartInstance = null;
    }
}

function drawSummaryChart(labels, data, colors) {
    destroySummaryChart();
    const canvas = document.getElementById("summaryChart");
    if (!canvas || typeof Chart === "undefined") return;
    const textColor = cssVar("--text", "#1e293b");
    const mutedColor = cssVar("--muted", "#64748b");
    const gridColor = cssVar("--border", "#e2e8f0");
    summaryChartInstance = new Chart(canvas, {
        type: "bar",
        data: {
            labels,
            datasets: [{
                data,
                backgroundColor: colors,
                borderColor: colors,
                borderWidth: 1,
                barThickness: "flex",
                maxBarThickness: 36,
            }],
        },
        options: {
            indexAxis: "y",
            responsive: true,
            maintainAspectRatio: false,
            layout: { padding: { right: 32, top: 8, bottom: 4, left: 4 } },
            plugins: {
                legend: { display: false },
                tooltip: {
                    callbacks: {
                        label: (ctx) => `${ctx.parsed.x} targets`,
                    },
                },
            },
            scales: {
                x: {
                    beginAtZero: true,
                    title: { display: true, text: "Number of targets", color: textColor },
                    ticks: { color: mutedColor, precision: 0 },
                    grid: { color: gridColor },
                },
                y: {
                    ticks: { color: textColor, autoSkip: false },
                    grid: { display: false },
                },
            },
        },
        plugins: [summaryValueLabelsPlugin],
    });
}

// Filter event wiring.
bystanderFilter.addEventListener("change", () => {
    if (currentRecords.length > 0) renderResultsTable();
});
scoreFilter.addEventListener("change", () => {
    if (currentRecords.length > 0) renderResultsTable();
});
if (gcFilter) {
    gcFilter.addEventListener("change", () => {
        if (currentRecords.length > 0) renderResultsTable();
    });
}
if (homopolymerFilter) {
    homopolymerFilter.addEventListener("change", () => {
        if (currentRecords.length > 0) renderResultsTable();
    });
}
if (trinucFilter) {
    trinucFilter.addEventListener("change", () => {
        if (currentRecords.length > 0) renderResultsTable();
    });
}
if (condenseToggle) {
    condenseToggle.addEventListener("change", () => {
        condenseEnabled = condenseToggle.checked;
        if (currentRecords.length > 0) renderResultsTable();
    });
}
