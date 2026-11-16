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
// Off-target screening (search off-targets for the filtered guide subset).
const screenOfftargetBtn = document.getElementById("screenOfftargetBtn");
const offtargetHint = document.getElementById("offtargetHint");
const otGenomeGroup = document.getElementById("otGenomeGroup");
const otGenomeSelect = document.getElementById("otGenomeSelect");
const otHideFailedToggle = document.getElementById("otHideFailedToggle");
const offtargetPanel = document.getElementById("offtargetPanel");
const offtargetPlaceholder = document.getElementById("offtargetPlaceholder");
const offtargetStatus = document.getElementById("offtargetStatus");
const downloadOfftargetResults = document.getElementById("downloadOfftargetResults");
const tabButtons = document.querySelectorAll(".tab-button");
const tabPanels = document.querySelectorAll(".tab-panel");
const subtabButtons = document.querySelectorAll(".subtab-button");
const subtabPanels = document.querySelectorAll(".subtab-panel");
const toastContainer = document.getElementById("toastContainer");
const modeSpreadsheetBtn = document.getElementById("modeSpreadsheet");
const modeSequenceBtn = document.getElementById("modeSequence");
const modeSingleBtn = document.getElementById("modeSingle");
const spreadsheetPanel = document.getElementById("spreadsheetPanel");
const sequencePanel = document.getElementById("sequencePanel");
const singlePanel = document.getElementById("singlePanel");
const sequenceInput = document.getElementById("sequenceInput");
const loadSequenceExampleBtn = document.getElementById("loadSequenceExampleBtn");
const loadSpreadsheetExampleBtn = document.getElementById("loadSpreadsheetExampleBtn");
const loadSingleExampleBtn = document.getElementById("loadSingleExampleBtn");
const genomeSelect = document.getElementById("genomeSelect");
const singleGenomeSelect = document.getElementById("singleGenomeSelect");
const singleChromosome = document.getElementById("singleChromosome");
const singlePosition = document.getElementById("singlePosition");
const singleEditType = document.getElementById("singleEditType");
const singleRunGroup = document.getElementById("singleRunGroup");
const defaultGenomeKey = genomeSelect ? genomeSelect.value : "human";
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
// The target set actually shown in the Input table and submitted to /run:
// the full upload when every target matches the selected editors, otherwise
// the filtered subset returned by /check (mismatching targets removed).
// Always re-derived from uploadedRawB64, which stays the untouched original
// so switching to a broader editor selection restores removed targets.
let activeRawB64 = null;
// Preview rows/columns of the ORIGINAL upload (before any compatibility
// filtering), kept so broadening the editor selection can restore the table.
let uploadedPreview = null;
let uploadedPreviewColumns = null;
let customGenomeB64 = null;
let currentJobId = null;
let currentRevokeToken = null;
let eventSource = null;
let pollHandle = null;
let referencePollHandle = null;
let currentRecords = [];
// When true, rows that differ only in `base_editor` are collapsed into a
// single row whose `base_editor` cell lists every editor newline-separated.
// Mirrors the `condense` query param on /download/results so the downloaded
// file matches what the user sees on screen. Initialized from the checkbox
// in the results toolbar (default: checked).
let condenseEnabled = !!(condenseToggle && condenseToggle.checked);
let matchedBaseEditorCount = 0;
// Off-target screening state. ``offtargetConfig`` mirrors the /config
// "offtarget" block (per-genome availability + operator thresholds);
// ``screenRecords``/``screenJobId`` hold the latest screening snapshot —
// unlike the Guides tab it is NOT re-filtered live, it shows exactly what
// was screened.
let offtargetConfig = null;
let screenJobId = null;
let screenRevokeToken = null;
let screenRecords = [];
let screening = false;
let currentFilteredGuideRows = [];
let currentJobGenomeKey = null;
// Input mode of the loaded job (not the currently selected input mode):
// drives the sequence-mode results display adaptation and the off-target
// search opt-in. Set when a run is submitted, alongside currentJobGenomeKey.
let currentJobMode = null;
// Whether the latest screening run was a sequence-mode search (whose
// on-target caveat is permanently shown below the search button); the
// off-target report repeats the caveat then.
let sequenceScreenAcknowledged = false;
// Verdicts of the latest off-target check: JSON screenRowKey ->
// {verdict: "pass"|"fail", reason}. Drives the Guides-table verdict column,
// the "hide failing guides" filter and the toggle title. Cleared together
// with the user's "hide failing guides" choice when a new pipeline run
// starts.
let otVerdictMap = new Map();
let otHideFailed = false;
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
        if (data.ok && data.offtarget) {
            offtargetConfig = data.offtarget;
            updateScreenButton();
            updateSearchInfo();
        }
    } catch (err) {
        console.warn("Failed to fetch server config:", err);
    }

    // Make sure the genome selector doesn't start on an unavailable option.
    ensureGenomeSelection();

    // Reference (genome + GTF) preload status: one-shot fetch, then poll
    // while loading.
    try {
        await refreshReferenceStatus();
    } catch (err) {
        console.warn("Failed to fetch reference status:", err);
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

// Combined reference (genome + GTF) preload status badge. Polls
// /reference/status while any preload is in flight, then stops on a
// terminal state (ready / error / disabled). A single aggregate pill is
// shown rather than one per genome; per-genome detail lives in the title.
async function refreshReferenceStatus() {
    const badge = document.getElementById("referenceStatus");
    if (!badge) return;

    let data;
    try {
        const res = await fetch("/reference/status");
        data = await res.json();
    } catch (err) {
        badge.textContent = "reference: status unavailable";
        badge.className = "badge error";
        return;
    }

    const detailLines = (data.genomes || []).map(g => {
        const parts = [`${g.key}: genome=${g.status}`];
        if (g.label) parts.push(`(${g.label})`);
        if (g.gtf_status !== undefined) {
            let gtf = `gtf=${g.gtf_status}`;
            if (g.gtf_label) gtf += ` (${g.gtf_label})`;
            parts.push(gtf);
        }
        return parts.join(" ");
    });
    badge.title = detailLines.join("\n");

    if (!data.ok || !data.preload_enabled) {
        badge.textContent = "reference: per-job";
        badge.className = "badge muted";
        stopReferencePoll();
        return;
    }

    if (data.status === "loading") {
        badge.textContent = "reference: loading\u2026";
        badge.className = "badge loading";
        if (!referencePollHandle) {
            referencePollHandle = setInterval(refreshReferenceStatus, 2000);
        }
        return;
    }

    if (data.status === "ready") {
        badge.textContent = "reference: ready";
        badge.className = "badge ready";
        stopReferencePoll();
        return;
    }

    if (data.status === "partial") {
        badge.textContent = "reference: partially loaded";
        badge.className = "badge error";
        stopReferencePoll();
        return;
    }

    if (data.status === "error") {
        badge.textContent = "reference: load failed";
        badge.className = "badge error";
        stopReferencePoll();
        return;
    }

    // not_loaded or anything unexpected — keep quiet but don't poll.
    badge.textContent = "reference: idle";
    badge.className = "badge muted";
    stopReferencePoll();
}

function stopReferencePoll() {
    if (referencePollHandle) {
        clearInterval(referencePollHandle);
        referencePollHandle = null;
    }
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

// Mode switcher: sequence, single and spreadsheet input. Each entry
// describes its tab button, input panel, where the shared Run button is
// anchored, whether Run starts enabled, and the Input-tab empty-state text.
const MODE_UI = {
    sequence: {
        button: modeSequenceBtn,
        panel: sequencePanel,
        runRow: () => sequencePanel.querySelector(".input-meta-row"),
        runEnabled: true,
        emptyState: "Paste sequences above and click Run to parse and start the pipeline.",
    },
    single: {
        button: modeSingleBtn,
        panel: singlePanel,
        runRow: () => singleRunGroup,
        runEnabled: true,
        emptyState: "Fill in the target fields above and click Run to parse and start the pipeline.",
    },
    spreadsheet: {
        button: modeSpreadsheetBtn,
        panel: spreadsheetPanel,
        runRow: () => document.getElementById("spreadsheetRunGroup"),
        runEnabled: false,
        emptyState: "Upload a targets file to preview input data.",
    },
};

function setMode(mode) {
    const active = MODE_UI[mode];
    if (!active) return;
    currentMode = mode;
    for (const [key, entry] of Object.entries(MODE_UI)) {
        const isActive = key === mode;
        entry.button.classList.toggle("active", isActive);
        entry.button.setAttribute("aria-selected", isActive ? "true" : "false");
        entry.button.tabIndex = isActive ? 0 : -1;
        entry.panel.style.display = isActive ? "" : "none";
    }

    // Move the shared Run button into the active input panel's action row so
    // it stays visually anchored next to the input controls. In spreadsheet
    // mode it lives in the right-aligned run group next to Load example; in
    // single mode it lives in the panel's run group next to the genome
    // selector; in sequence mode it sits next to the Load example button,
    // whose 'right' class pushes the pair to the far edge of the row.
    const targetRow = active.runRow();
    if (targetRow) {
        targetRow.appendChild(runBtn);
    }

    // Reset state when switching modes: a previously uploaded file is
    // irrelevant for other inputs and vice-versa. In sequence and single
    // mode the Run button starts enabled because clicking it triggers the
    // parse + check + run flow itself; in spreadsheet mode it stays
    // disabled until a file is uploaded (handled by uploadTargetsFile →
    // runCompatibilityCheck).
    uploadedRawB64 = null;
    activeRawB64 = null;
    uploadedPreview = null;
    uploadedPreviewColumns = null;
    customGenomeB64 = null;
    matchedBaseEditorCount = 0;
    runBtn.disabled = !active.runEnabled;
    compatibilityInfo.textContent = "";
    checkAlert.style.display = "none";
    inputAlert.style.display = "none";
    inputPanel.innerHTML = `<p style="color: var(--muted);padding:16px;">${escapeHtml(active.emptyState)}</p>`;
}

modeSpreadsheetBtn.addEventListener("click", () => setMode("spreadsheet"));
modeSequenceBtn.addEventListener("click", () => setMode("sequence"));
modeSingleBtn.addEventListener("click", () => setMode("single"));

// Pressing Enter in the single-target fields starts the run flow.
for (const field of [singleChromosome, singlePosition]) {
    if (field) {
        field.addEventListener("keydown", (ev) => {
            if (ev.key === "Enter" && !runBtn.disabled) {
                ev.preventDefault();
                runBtn.click();
            }
        });
    }
}

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
    // Always check the ORIGINAL upload, never a previously filtered subset,
    // so switching to a broader editor selection brings removed targets back.
    form.append("raw_b64", uploadedRawB64);
    form.append("base_editor", baseEditorSelectionValue());
    try {
        const res = await fetch("/check", { method: "POST", body: form });
        const data = await res.json();
        if (!data.ok) return;

        if (data.matched === 0) {
            // No target can be processed: keep the hard error and the full
            // preview table visible, and refuse to run.
            showAlert(checkAlert, "error",
                `${data.mismatched}/${data.total} target(s) do not match the selected base editor profile(s). ` +
                `Supported types: ${data.editor_types.join(", ")}. ` +
                `Please choose a different profile or fix your input.`
            );
            compatibilityInfo.textContent = `0/${data.total} targets match the selected base editor profile(s)`;
            compatibilityInfo.style.color = "var(--danger)";
            matchedBaseEditorCount = 0;
            activeRawB64 = null;
            if (currentMode !== "sequence" && uploadedPreview) {
                // Restore the full table in case a previous mixed check had
                // narrowed it to the acceptable subset.
                inputPanel.innerHTML = renderTable(uploadedPreview, uploadedPreviewColumns);
            }
            runBtn.disabled = true;
            return;
        }

        // At least one target matches. When some targets mismatch, replace
        // the input table and the /run payload with the acceptable subset.
        // In sequence mode the input panel shows a per-sequence preview card
        // rather than the target rows, so leave it untouched there.
        activeRawB64 = data.mismatched > 0 ? data.filtered_raw_b64 : uploadedRawB64;
        matchedBaseEditorCount = data.matched;
        if (data.mismatched > 0) {
            if (currentMode !== "sequence") {
                inputPanel.innerHTML = renderTable(data.filtered_preview, data.filtered_columns);
            }
            showAlert(checkAlert, "warning",
                `${data.mismatched} of ${data.total} target(s) ` +
                (data.mismatched === 1 ? "does" : "do") +
                ` not have acceptable base changes for the selected base editor ` +
                `profile(s) (supported: ${data.editor_types.join(", ")}) ` +
                `and ${data.mismatched === 1 ? "was" : "were"} removed from the input. ` +
                `The remaining ${data.matched} target(s) will be used for the run.`
            );
        } else {
            if (currentMode !== "sequence" && uploadedPreview) {
                // Everything matches again (e.g. the selection was broadened):
                // restore the original, unfiltered table.
                inputPanel.innerHTML = renderTable(uploadedPreview, uploadedPreviewColumns);
            }
            checkAlert.style.display = "none";
        }
        compatibilityInfo.textContent = `${data.matched}/${data.total} targets match the selected base editor profile(s)`;
        compatibilityInfo.style.color = "var(--success)";
        runBtn.disabled = false;
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

// Single-mode parse helper. Posts the single-target form fields to
// /single/parse and returns a structured result. No UI side effects —
// callers decide how to surface errors and whether to switch tabs.
async function parseSingleTarget() {
    const chromosome = singleChromosome ? singleChromosome.value.trim() : "";
    const position = singlePosition ? singlePosition.value.trim() : "";
    const edit = singleEditType ? singleEditType.value : "";
    if (!chromosome || !position || !edit) {
        return { ok: false, error: "Fill in chromosome, position and editing type." };
    }
    const form = new FormData();
    form.append("chromosome", chromosome);
    form.append("position", position);
    form.append("edit", edit);
    form.append("genome", selectedGenomeKey());
    let res, data;
    try {
        res = await fetch("/single/parse", { method: "POST", body: form });
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
        preview: data.preview,
        columns: data.columns,
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

// Return the genome key currently selected for the active input mode. The
// spreadsheet and single panels each have their own genome selector; the
// default is used when the selector is absent.
function selectedGenomeKey() {
    const select = currentMode === "single" ? singleGenomeSelect : genomeSelect;
    if (select && select.value) return select.value;
    return defaultGenomeKey || "human";
}

// If a server-rendered default genome is unavailable (its <option> is
// disabled because the FASTA is missing on disk), fall back to the first
// available genome so no selector ever submits a reference the server
// cannot load.
function ensureGenomeSelection() {
    for (const select of [genomeSelect, singleGenomeSelect]) {
        if (!select) continue;
        const current = select.selectedOptions && select.selectedOptions[0];
        if (current && !current.disabled) continue;
        const fallback = Array.from(select.options).find(o => !o.disabled);
        if (fallback) select.value = fallback.value;
    }
}

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
    form.append("genome", selectedGenomeKey());

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
        activeRawB64 = null;
        uploadedPreview = null;
        uploadedPreviewColumns = null;
        matchedBaseEditorCount = 0;
        runBtn.disabled = true;
        compatibilityInfo.textContent = "";
        checkAlert.style.display = "none";
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(data.error || "Upload failed.")}</p>`;
        return;
    }

    showAlert(inputAlert, "success", `Loaded ${data.row_count} rows from "${file.name}". Showing up to 100.`);
    uploadedRawB64 = data.raw_b64;
    // Optimistic: show the full upload; runCompatibilityCheck (below) narrows
    // the table to the acceptable subset when some targets mismatch.
    activeRawB64 = data.raw_b64;
    uploadedPreview = data.preview;
    uploadedPreviewColumns = data.columns;
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
        // The example targets are human (GRCh38); force the selector to
        // human so chromosome validation uses the right reference even when
        // the operator's default genome is mouse.
        if (genomeSelect) {
            const humanOption = Array.from(genomeSelect.options)
                .find(o => o.value === "human");
            if (humanOption) genomeSelect.value = "human";
        }
        const file = new File(
            [SPREADSHEET_EXAMPLE_TSV],
            "example_targets.tsv",
            { type: "text/tab-separated-values" }
        );
        await uploadTargetsFile(file);
    });
}

// Single-mode "Load example" button. Fills the fields with the first target
// of the spreadsheet example (human GRCh38) so the example is consistent
// across modes. Like the sequence-mode example it only fills the input —
// the run itself stays a separate Run click.
const SINGLE_EXAMPLE = { chromosome: "12", position: "25245350", edit: "T>C" };

if (loadSingleExampleBtn) {
    loadSingleExampleBtn.addEventListener("click", () => {
        if (currentMode !== "single") setMode("single");
        singleChromosome.value = SINGLE_EXAMPLE.chromosome;
        singlePosition.value = SINGLE_EXAMPLE.position;
        singleEditType.value = SINGLE_EXAMPLE.edit;
        // The example coordinates are human (GRCh38); force the panel's
        // selector to human so validation uses the right reference even
        // when the operator's default genome is mouse.
        if (singleGenomeSelect) {
            const humanOption = Array.from(singleGenomeSelect.options)
                .find(o => o.value === "human");
            if (humanOption) singleGenomeSelect.value = "human";
        }
    });
}

// Run pipeline
runBtn.addEventListener("click", async () => {
    if (selectedEditors.size === 0) {
        showToast("error", "Select at least one base editor.");
        return;
    }

    // Sequence mode: parse the textarea on every Run click so the run always
    // uses the current input (re-parsing avoids stale submissions after
    // edits). On parse error, route the user to the Input tab so the error
    // is visible; on success, render the preview silently into inputPanel
    // and continue without switching tabs.
    let parsedThisClick = false;
    if (currentMode === "sequence") {
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
        // Fresh input invalidates any previously derived (possibly filtered)
        // subset; runCompatibilityCheck below recomputes it. Resetting here
        // also guarantees a failed check aborts the run instead of letting a
        // stale subset from a previous run slip through.
        activeRawB64 = null;
        inputPanel.innerHTML = renderSequencePreview(parsed.preview);
        parsedThisClick = true;
    }

    // Single mode: parse the form fields on every Run click (cheap for one
    // row, and re-parsing avoids stale submissions after edits). Error
    // handling mirrors the sequence path; success renders the standard
    // 1-row target table into inputPanel.
    if (currentMode === "single") {
        inputAlert.style.display = "none";
        checkAlert.style.display = "none";
        compatibilityInfo.textContent = "";
        inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">Parsing target\u2026</p>`;
        const parsed = await parseSingleTarget();
        if (!parsed.ok) {
            switchTab("inputTab");
            showAlert(inputAlert, "error", parsed.error);
            inputPanel.innerHTML = `<p style="color:var(--muted);padding:16px;">${escapeHtml(parsed.error)}</p>`;
            runBtn.disabled = false;
            return;
        }
        uploadedRawB64 = parsed.raw_b64;
        activeRawB64 = parsed.raw_b64;
        uploadedPreview = parsed.preview;
        uploadedPreviewColumns = parsed.columns;
        inputPanel.innerHTML = renderTable(parsed.preview, parsed.columns);
        parsedThisClick = true;
    }

    if (!uploadedRawB64) return;

    // When parse just happened, re-run the compatibility check before
    // submitting so the user is alerted to any target/editor mismatches. On
    // mismatch, route to the Input tab where the checkAlert is shown. The
    // check also derives activeRawB64 (the acceptable subset) which is what
    // gets submitted below.
    if (parsedThisClick) {
        await runCompatibilityCheck();
        if (runBtn.disabled) {
            switchTab("inputTab");
            return;
        }
    }
    if (!activeRawB64) return;

    runBtn.disabled = true;
    setStatus("submitting");
    clearLogs();
    currentRecords = [];
    resultsPanel.querySelectorAll(".results-table-wrap").forEach(n => n.remove());
    if (resultsPlaceholder) {
        resultsPlaceholder.style.display = "block";
        // Reset the colour in case a previous run's "No results returned."
        // message left the placeholder in the dark red "nothing found" tone.
        resultsPlaceholder.style.color = "var(--muted)";
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
    // Any previous off-target snapshot belongs to an older parent job.
    resetScreeningState();

    // Switch to Logging tab
    tabButtons.forEach(b => b.classList.remove("active"));
    tabPanels.forEach(p => p.classList.remove("active"));
    document.querySelector('[data-tab="logTab"]').classList.add("active");
    document.getElementById("logTab").classList.add("active");

    const form = new FormData();
    form.append("raw_b64", activeRawB64);
    form.append("base_editor", baseEditorSelectionValue());
    form.append("mode", currentMode);
    form.append("genome", selectedGenomeKey());
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
    // Remember which genome this job runs against (so the off-target
    // screening button can check availability for the right index) and
    // which input mode it came from (drives the sequence-mode results
    // display and the screening opt-in, even if the user switches input
    // modes while the results are open).
    currentJobGenomeKey = selectedGenomeKey();
    currentJobMode = currentMode;
    updateSearchInfo();
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

// Off-target screening: search off-targets for the currently filtered
// subset of guides (button in the Results-tab filter bar).
if (screenOfftargetBtn) {
    screenOfftargetBtn.addEventListener("click", runScreening);
}

// Sequence-mode reference-genome selector: availability (and the info
// tooltip on the section title) depend on the chosen genome.
if (otGenomeSelect) {
    otGenomeSelect.addEventListener("change", () => {
        updateScreenButton();
        updateSearchInfo();
    });
}

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
// value renders as an em-dash. ``homopolymer`` is rendered as plain text
// ("poly" when flagged, "ok" otherwise) to match the rest of the table.
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
            return "poly";
        }
        if (value === false || value === "false") {
            return "ok";
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

// Render the off-target summary columns produced by the BEguidedOFF step.
// Guides whose worst reported off-target score is empty had no off-target at
// or below the search cutoff — they pass the off-target filter — so the cell
// says so instead of rendering a blank (the underlying TSV keeps an empty
// cell; only the on-screen table shows the label).
//
// When operator thresholds are known (from /config), values that the
// pipeline's removal rule would have removed (score below the threshold or
// too many off-target sites) get an explicit badge. The thresholds are
// display-only: users never configure them.
// Evaluate one screened guide row against the operator's off-target
// thresholds (the same rule the strict pipeline filter applies). Returns
// null when the row carries no screening result, otherwise
// {verdict: "pass"|"fail", reason: "<human-readable why>"}.
function offtargetVerdict(row) {
    const thresholds = offtargetConfig && offtargetConfig.thresholds;
    if (!thresholds) return null;
    const has = (v) => v !== undefined && v !== null && v !== "";
    const scoreRaw = row.worst_offtarget_score;
    const countRaw = row.offtarget_count;
    if (!has(scoreRaw) && !has(countRaw)) return null;
    if (has(scoreRaw)) {
        const num = Number(scoreRaw);
        if (Number.isFinite(num) && num < thresholds.remove_guides_score) {
            return {
                verdict: "fail",
                reason: `severe off-target: worst score ${num.toFixed(2)} is below ` +
                    `the operator threshold (${thresholds.remove_guides_score})`,
            };
        }
    }
    if (has(countRaw)) {
        const num = Number(countRaw);
        if (Number.isFinite(num) && num > thresholds.remove_guides_n_offtargets) {
            return {
                verdict: "fail",
                reason: `too many off-targets: ${num} sites exceed the ` +
                    `operator threshold (${thresholds.remove_guides_n_offtargets})`,
            };
        }
    }
    return {
        verdict: "pass",
        reason: "no severe off-target within the search limit; " +
            "milder off-targets may exist",
    };
}

// Pass/fail badge with the reason as tooltip; shared by the Guides-table
// verdict column and the compact report tab.
function verdictBadgeHtml(v) {
    if (!v) return `<span class="muted">&mdash;</span>`;
    if (v.verdict === "fail") {
        return `<span class="ot-badge ot-badge-danger" title="${escapeHtml(v.reason)}">fail</span>`;
    }
    return `<span class="ot-badge ot-badge-pass" title="${escapeHtml(v.reason)}">pass</span>`;
}

// Verdict cell for the compact report tab: one row per screened guide, so
// the verdict comes straight from the row's own screening columns.
function formatReportVerdictCell(row) {
    return verdictBadgeHtml(offtargetVerdict(row));
}

// Verdict cell for the Guides table. Rows carry an attached
// `offtarget_verdict` value ("" = not part of the screened subset); the
// tooltip lists the per-editor detail, which matters for condensed rows
// whose editors share one verdict.
function formatVerdictCell(row) {
    if (!screenRecords.length) return `<span class="muted">&mdash;</span>`;
    const v = row.offtarget_verdict;
    if (v === undefined || v === null || v === "") {
        if (isNoGuideRow(row)) return `<span class="muted">&mdash;</span>`;
        return `<span class="muted" title="This guide was not part of the last off-target check.">not screened</span>`;
    }
    const editors = String(row.base_editor ?? "")
        .split("\n").map(s => s.trim()).filter(Boolean);
    const parts = [];
    for (const ed of editors) {
        const res = otVerdictMap.get(JSON.stringify(screenRowKey({ ...row, base_editor: ed })));
        parts.push(`${ed || "(editor)"}: ` +
            (res ? `${res.verdict} \u2014 ${res.reason}` : "not screened"));
    }
    const tip = parts.join("\n") || v;
    const badge = v === "fail"
        ? `<span class="ot-badge ot-badge-danger" title="${escapeHtml(tip)}">fail</span>`
        : `<span class="ot-badge ot-badge-pass" title="${escapeHtml(tip)}">pass</span>`;
    return badge;
}

function formatOfftargetCell(row, column, value) {
    const thresholds = offtargetConfig && offtargetConfig.thresholds;
    if (column === "worst_offtarget_score") {
        if (value === undefined || value === null || value === "") {
            return `<span class="muted">none within search limit</span>`;
        }
        const num = Number(value);
        if (Number.isFinite(num)) {
            let html = num.toFixed(2).replace(/\.?0+$/, "");
            if (thresholds && num < thresholds.remove_guides_score) {
                html += ` <span class="ot-badge ot-badge-danger" title="Worst off-target score is below the operator's off-target threshold (${thresholds.remove_guides_score}): this guide would be removed by the strict off-target filter. Lower CCTOP score = more severe.">severe off-target</span>`;
            }
            return html;
        }
        return null;
    }
    if (column === "offtarget_count") {
        if (value === undefined || value === null || value === "") {
            return `<span class="muted">none within search limit</span>`;
        }
        const num = Number(value);
        if (Number.isFinite(num) && thresholds &&
            num > thresholds.remove_guides_n_offtargets) {
            return `${num} <span class="ot-badge ot-badge-danger" title="More off-target sites than the operator's threshold (${thresholds.remove_guides_n_offtargets}); this guide would be removed by the strict off-target filter.">too many</span>`;
        }
        return null;
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
    offtarget_verdict: "Off-target check",
    trinuc_edit_score: "context edit score",
    worst_offtarget_score: "worst off-target score",
    offtarget_count: "off-target count",
};

// Tooltip text shown on hover for results-table headers. Keyed by the
// snake_case column name (as emitted by the pipeline), so the lookup works
// regardless of the display label. Only the results-table call site passes
// this map to renderTable(); the upload-preview table is left untouched.
const RESULTS_COLUMN_TOOLTIPS = {
    gene_edit_desired:
        "Desired edit at the protein level, formatted as '<gene> <Target_AA><codon_number><Desired_AA>' (e.g. 'TP53 *306R'). Codon numbering follows the MANE Select CDS frame; empty when the target is not in a CDS or no gene name resolves.",
    genome_base:
        "Reference-genome base at the target position. Useful for confirming that the target position and edit are correct.",
    editing_position:
        "1-based position of the edited base, counted from the 5' end of the guide RNA (not the genomic coordinate).",
    bystander_positions:
        "Comma-separated 1-based positions from the 5' end of the guide where additional bases may be edited as bystanders.",
    bystander_effect:
        "Protein-level consequence of the bystander edits: synonymous (same amino acid), nonsynonymous (any amino-acid change, including missense and nonsense), or none (no bystander edits). Assessed against the MANE Select transcript codon frame. Empty when no GTF was provided.",
    guide_score:
        "CRISPRater on-target activity score (0-1). Derived from spCas9 training data and used here as a general guide-quality proxy.",
    gc_content:
        "GC fraction of the guide RNA spacer, expressed as a percentage. Guides with very low or very high GC tend to bind less efficiently.",
    homopolymer:
        "True when the spacer contains a run of >=4 identical bases (e.g. AAAA). Homopolymer runs can hurt guide specificity and efficiency.",
    trinuc_edit_score:
        "Marquart et al. 2021 normalised ABE8e editing proportion for the local X-A-Y trinucleotide context (the bases immediately before and after the target A). 0-1; NaN for CBE guides.",
    base_editor:
        "Generally, favor SpCas9 over NG/SpG over SpRY. Relaxed PAMs trade activity and specificity for targeting range. Regarding ABEs, ABE8e is most established and a reliable choice.",
    offtarget_verdict:
        "Verdict of the off-target check against the operator's thresholds. The search enumerates only the most severe off-targets (CCTOP score at or below the operator's score threshold), so 'pass' means none that severe were found \u2014 guides can still have milder off-targets. Hover a badge for the specific reason.",
    worst_offtarget_score:
        "Lowest CCTOP score among the guide's off-targets; lower means more concerning (a perfect-match site scores 0.0). Empty means no off-target at or below the search cutoff was found. Guides whose worst score falls below the removal threshold are removed by the strict off-target filter.",
    offtarget_count:
        "Number of editable off-targets with a CCTOP score at or below the count threshold (fast mode counts at the search cutoff). Guides with more than the allowed maximum are removed from the results.",
    worst_offtarget_mismatch_pos:
        "Comma-separated 1-based positions (from the 5' end of the guide) where the worst off-target differs from the guide RNA.",
    worst_offtarget_chr:
        "Chromosome of the worst (most concerning) off-target for this guide.",
    worst_offtarget_pos:
        "1-based genomic position of the worst off-target (leftmost base for + strand guides, rightmost for - strand guides).",
    worst_offtarget_strand:
        "Strand of the worst off-target: '+' if it matches the guide's genomic strand, '-' for the reverse complement.",
    worst_offtarget_sequence:
        "Genomic sequence (20 bp spacer) of the worst off-target, oriented to align with the guide RNA.",
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
    if (currentJobMode === "sequence"
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
            // "Nothing found" states stand out in dark red instead of grey.
            resultsPlaceholder.style.color = "#b91c1c";
            resultsPlaceholder.textContent = data.message || "No results returned.";
        }
        resultsToolbar.style.display = "none";
        if (resultsSummary) {
            resultsSummary.innerHTML = "";
            resultsSummary.style.display = "none";
        }
        setSummaryPlaceholder(true, data.message || "No results returned.", "danger");
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
        if (mode === "none") {
            if (effect) return effect === "none";
            const positions = (row.bystander_positions ?? "").toString().trim();
            return positions === "";
        }
        if (mode === "synonymous") {
            if (effect === "synonymous" || effect === "none") return true;
            // No bystander metadata (e.g. sequence mode, where step 2 is
            // skipped and the column is absent): treat the guide as
            // bystander-free, mirroring the "none" branch above.
            if (!effect) {
                const positions = (row.bystander_positions ?? "").toString().trim();
                return positions === "";
            }
            return false;
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
// result column is emptied. Records keep the spreadsheet identity fields
// (chromosome/position, target_base/desired_base) in every mode;
// ``sequence`` and ``desired_edit`` are legacy display aliases kept as
// harmless fallbacks. Mirrors filters._PLACEHOLDER_IDENTITY_COLUMNS.
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

// After filtering, order the rows by input order and interleave one
// placeholder row per target (unique ``targetKey``) that has no surviving
// row. The pre-filter records carry the pipeline's input order, so walking
// them puts every target back at its original position — surviving rows keep
// their current order within the target. Mirrors
// ``beguided.filtering.filters._fill_missing_target_rows`` so the on-screen
// table matches the downloaded file.
function fillMissingTargets(filtered, preFilter) {
    if (!preFilter || preFilter.length === 0) return filtered;
    // Group the surviving rows by target, preserving their current order.
    const grouped = new Map();
    for (const row of filtered) {
        const key = targetKey(row);
        if (!grouped.has(key)) grouped.set(key, []);
        grouped.get(key).push(row);
    }
    const out = [];
    const seen = new Set();
    for (const row of preFilter) {
        const key = targetKey(row);
        if (seen.has(key)) continue;
        seen.add(key);
        const rows = grouped.get(key);
        if (rows && rows.length) {
            out.push(...rows);
        } else {
            // First pre-filter record of the target is the template for the
            // placeholder (it carries the identity/base info even when the
            // target's guides were all removed).
            out.push(makeNoGuideRow(row));
        }
    }
    // Safety net: never drop a surviving row whose target is absent from the
    // pre-filter set (cannot normally happen).
    for (const [key, rows] of grouped) {
        if (!seen.has(key)) out.push(...rows);
    }
    return out;
}

// Attach the screening verdict to every results row before rendering and
// condensing. Returns new row objects (the shared currentRecords stay
// untouched); "" marks guides that were not part of the screened subset.
// The verdict is attached BEFORE condensing so it becomes part of the
// condense group key: a target whose editors disagree then splits into a
// failing row and a passing row, each listing only same-verdict editors.
// With no snapshot the records pass through unchanged and the Guides table
// shows no verdict column at all.
function attachVerdicts(records) {
    if (!screenRecords.length || !otVerdictMap.size) return records;
    return records.map(r => {
        const v = otVerdictMap.get(JSON.stringify(screenRowKey(r)));
        return { offtarget_verdict: v ? v.verdict : "", ...r };
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

    const { mode, score, gc, homopolymer, trinuc } = getActiveFilters();
    // Attach the screening verdict (if a snapshot exists) before condensing
    // so same-verdict rows merge and mixed targets split into separate
    // pass/fail rows; then condense BEFORE filtering so the filter sees the
    // same row shape the download endpoint will return. When there is no
    // snapshot and the toggle is off both steps pass records through
    // unchanged.
    const annotated = attachVerdicts(currentRecords);
    const condensed = condenseEnabled ? condenseRecords(annotated) : annotated;
    const attributeFiltered = filterRecords(condensed, mode, score, gc, homopolymer, trinuc);
    // Re-add one placeholder row per target whose guides were all removed by
    // the active filters (mirrors the server-side download fill).
    const placeholderFilled = fillMissingTargets(attributeFiltered, currentRecords);
    // Optional client-side "Hide failing guides" filter: removes rows whose
    // verdict is fail. One rule serves both table modes because condensed
    // rows are verdict-homogeneous — hiding a failing condensed row removes
    // only the failing editors while the target stays visible through its
    // passing sub-row. Applied AFTER the placeholder fill (so it never
    // re-adds misleading "**No guide possible**" placeholders for targets
    // whose guides all failed — those guides exist, they just failed) and
    // never hides the placeholder rows themselves. Nothing is deleted:
    // untick the box to see every guide again.
    const hideFailed = otHideFailed && otVerdictMap.size > 0;
    const keepRow = (r) => isNoGuideRow(r) || r.offtarget_verdict !== "fail";
    const filtered = hideFailed ? placeholderFilled.filter(keepRow) : placeholderFilled;
    // Summary is intentionally NOT affected by the condense toggle: compute a
    // separate filtered view from the uncondensed records so per-editor counts
    // always reflect one editor per row. Filters still apply, including the
    // off-target "hide failing guides" filter.
    const summaryUnfiltered = filterRecords(annotated, mode, score, gc, homopolymer, trinuc);
    const summaryFiltered = hideFailed ? summaryUnfiltered.filter(keepRow) : summaryUnfiltered;

    // Off-target screening counts the real guide rows (not the re-added
    // "no guide" placeholders) the current filters select — including the
    // hide-failing filter, so the button count matches the visible guides.
    currentFilteredGuideRows = filtered.filter(r => !isNoGuideRow(r));
    updateScreenButton();

    target.classList.toggle("condensed", condenseEnabled && condensed.length < currentRecords.length);

    if (filtered.length === 0) {
        target.innerHTML = `<p style="color:#b91c1c;padding:16px;">No guides match the current filters.</p>`;
    } else {
        const columns = Object.keys(filtered[0]);
        // With a screening snapshot every annotated row carries the verdict,
        // but re-added placeholder rows may lack the key; keep the column
        // when the first visible row happens to be a placeholder.
        if (screenRecords.length && !columns.includes("offtarget_verdict")) {
            columns.unshift("offtarget_verdict");
        }
        // Sequence-mode display adaptation. The records keep the raw
        // screening-relevant fields (chromosome, position, target_base,
        // desired_base) so off-target screening keys match, but the table
        // mirrors the input-oriented view: the sequence name first, the
        // edit shown as a single "desired edit" column, and coordinates
        // that only exist relative to the pasted sequence hidden.
        let finalColumns;
        if (currentJobMode === "sequence") {
            let cols = columns.filter(
                c => c !== "position" && c !== "strand_used" && c !== "genome_base"
            );
            if (cols.includes("target_base")) {
                const at = cols.indexOf("target_base");
                cols = cols.filter(c => c !== "target_base" && c !== "desired_base");
                cols.splice(at, 0, "desired_edit");
            }
            finalColumns = reorderColumns(cols, "chromosome");
        } else {
            finalColumns = columns;
        }
        const displayColumns = finalColumns.map(columnLabel);
        target.innerHTML = renderTable(
            filtered,
            finalColumns,
            displayColumns,
            targetKey,
            (row, column, value) => {
                if (column === "offtarget_verdict") return formatVerdictCell(row);
                if (column === "desired_edit") {
                    return escapeHtml(`${row.target_base ?? ""}>${row.desired_base ?? ""}`);
                }
                const seqProps = formatSequencePropsCell(row, column, value);
                if (seqProps !== null) return seqProps;
                const otCell = formatOfftargetCell(row, column, value);
                if (otCell !== null) return otCell;
                return formatDnaSequenceCell(row, column, value);
            },
            resultsColumnTooltips(),
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

// ----- Off-target screening (filtered subset) ------------------------------
//
// Two-phase model: the guide search runs the pipeline without BEguidedOFF;
// once the user has narrowed the Guides table with the filters above, the
// "Check severe off-targets" button screens exactly the visible guide rows.
// The verdict lands as an extra first column in the Guides table itself
// (joined client-side from the immutable snapshot), and the "Off-target
// report" sub-tab holds a compact per-guide specificity report with the
// worst-site details.

// Identity of one results row for screening requests. Must match the
// server-side ``_screen_row_key`` in beguided/web/tasks.py: chromosome,
// position, target_base, desired_base and base_editor together uniquely
// identify a guide row of a completed job.
function screenRowKey(row) {
    const norm = (v) => (v === undefined || v === null) ? "" : String(v);
    return {
        chromosome: norm(row.chromosome),
        position: norm(row.position),
        target_base: norm(row.target_base),
        desired_base: norm(row.desired_base),
        base_editor: norm(row.base_editor),
    };
}

// Expand the currently filtered guide rows into per-guide screening keys.
// Condensed rows (one row per target, editors newline-joined) expand to one
// key per editor, so every visible guide is screened exactly once.
function collectScreeningRows(rows) {
    const out = [];
    const seen = new Set();
    for (const row of rows || []) {
        if (isNoGuideRow(row)) continue;
        const editors = String(row.base_editor ?? "")
            .split("\n").map(s => s.trim()).filter(Boolean);
        const editorList = editors.length ? editors : [""];
        for (const editor of editorList) {
            const key = screenRowKey({ ...row, base_editor: editor });
            const k = JSON.stringify(key);
            if (seen.has(k)) continue;
            seen.add(k);
            out.push(key);
        }
    }
    return out;
}

// Reset the screening snapshot (called when a new pipeline run starts — the
// old snapshot belongs to a different parent job).
function resetScreeningState() {
    screenJobId = null;
    screenRevokeToken = null;
    screenRecords = [];
    screening = false;
    sequenceScreenAcknowledged = false;
    currentFilteredGuideRows = [];
    otVerdictMap = new Map();
    otHideFailed = false;
    if (otHideFailedToggle) {
        otHideFailedToggle.checked = false;
        otHideFailedToggle.disabled = true;
        otHideFailedToggle.title = "";
    }
    const checkLine = document.getElementById("otCheckStatus");
    if (checkLine) checkLine.remove();
    if (offtargetPanel) {
        offtargetPanel.querySelectorAll(".results-table-wrap").forEach(n => n.remove());
    }
    if (offtargetPlaceholder) {
        offtargetPlaceholder.style.display = "block";
        offtargetPlaceholder.textContent =
            "The off-target report appears here after you run the off-target check on the filtered guides (button in the filter bar above the results).";
    }
    if (offtargetStatus) {
        offtargetStatus.style.display = "none";
        offtargetStatus.textContent = "";
    }
    if (downloadOfftargetResults) downloadOfftargetResults.style.display = "none";
    updateScreenButton();
}

// Human-readable description of the fast-mode search window: the exact
// PAM-proximal seed and the tolerated PAM-distal mismatches. The numbers
// come from /config (derived from the same cutoff the screening uses) so
// the text always matches the engine.
function searchWindowText(fs) {
    return (
        `The search enumerates only off-targets that match the guide exactly ` +
        `over the ${fs.seed_bp} guide bases closest to the PAM (positions ` +
        `${fs.seed_span}); mismatches are tolerated only in the PAM-distal ` +
        `${fs.hot_bp} bases (up to ${fs.hot_mismatches}). Sites with a ` +
        `mismatch inside that region score above the search cutoff and are ` +
        `not enumerated.`
    );
}

// Verdict-column header tooltip: what the verdict means plus the search
// window behind it. Falls back to the static generic text when the fast
// search summary is unavailable.
function verdictColumnTooltip() {
    const fs = offtargetConfig && offtargetConfig.fast_search;
    if (!fs) return RESULTS_COLUMN_TOOLTIPS.offtarget_verdict;
    return (
        `Verdict of the off-target check against the operator's thresholds. ` +
        searchWindowText(fs) +
        ` 'pass' therefore means no off-target that severe was found \u2014 ` +
        `guides can still have milder off-targets. Hover a badge for the ` +
        `specific reason.`
    );
}

// Column header tooltips for the results tables: the static map with the
// verdict-column tooltip refreshed from the live search-window facts.
function resultsColumnTooltips() {
    const tips = { ...RESULTS_COLUMN_TOOLTIPS };
    tips.offtarget_verdict = verdictColumnTooltip();
    return tips;
}

// Info tooltip on the "Off-target check" section title: the same search
// window description (native title attribute on the title, info glyph as
// the visual affordance — same pattern as the filter labels). The glyph is
// hidden when the fast engine is absent rather than showing a possibly
// wrong number.
function updateSearchInfo() {
    const titleEl = document.getElementById("otSearchTitle");
    if (!titleEl) return;
    const glyph = titleEl.querySelector(".th-info");
    const fs = offtargetConfig && offtargetConfig.fast_search;
    const genomes = (offtargetConfig && offtargetConfig.genomes) || {};
    const key = screeningGenomeKey();
    const mode = genomes[key] && genomes[key].mode;
    if (!fs || mode !== "fast") {
        titleEl.removeAttribute("title");
        if (glyph) glyph.style.display = "none";
        return;
    }
    titleEl.title = searchWindowText(fs);
    if (glyph) glyph.style.display = "";
}

// The genome the off-target search will run against: in sequence mode the
// user picks it in the selector at the right of the Off-target check
// controls; otherwise it is the genome the job ran with.
function screeningGenomeKey() {
    if (currentJobMode === "sequence" && otGenomeSelect && otGenomeSelect.value) {
        return otGenomeSelect.value;
    }
    return currentJobGenomeKey || defaultGenomeKey;
}

// Sequence-mode caveat, shown below the search button: without genomic
// coordinates the on-target site cannot be excluded from the search.
const SEQUENCE_OFFTARGET_NOTE =
    "Sequence mode: the on-target site cannot be excluded because the " +
    "input sequence's location in the reference genome is unknown. If " +
    "your sequence comes from the searched genome, the worst hit may be " +
    "the on-target itself (a perfect match). For plasmids or inserts " +
    "that are not in the genome, all hits are genuine off-targets.";

// Why the screening button cannot be used right now (null = it can).
function screeningBlockReason(nGuides) {
    if (!offtargetConfig) return "server configuration unavailable";
    const thresholdsAvail = offtargetConfig.genomes || {};
    const key = screeningGenomeKey();
    const avail = thresholdsAvail[key];
    if (avail && !avail.available) {
        return avail.reason || "off-target search not configured for this genome";
    }
    const cap = offtargetConfig.max_screen_guides;
    if (nGuides === 0) return "no guides match the current filters";
    if (cap && nGuides > cap) {
        return `too many guides selected (${nGuides}); narrow the filters below ${cap} and retry`;
    }
    return null;
}

// Reflect N = number of filtered guides + all disabled reasons on the button.
function updateScreenButton() {
    if (!screenOfftargetBtn) return;
    // The reference-genome selector only applies to sequence jobs (the
    // other modes inherit the genome chosen in the input panel).
    if (otGenomeGroup) {
        otGenomeGroup.style.display = currentJobMode === "sequence" ? "" : "none";
    }
    if (screening) {
        screenOfftargetBtn.disabled = true;
        screenOfftargetBtn.classList.add("btn-loading");
        screenOfftargetBtn.setAttribute("aria-busy", "true");
        screenOfftargetBtn.textContent = "Screening\u2026";
        screenOfftargetBtn.title = "Off-target check is running.";
        if (offtargetHint) {
            offtargetHint.textContent = "";
            offtargetHint.style.display = "none";
        }
        return;
    }
    screenOfftargetBtn.classList.remove("btn-loading");
    screenOfftargetBtn.removeAttribute("aria-busy");
    const rows = collectScreeningRows(currentFilteredGuideRows);
    const reason = screeningBlockReason(rows.length);
    screenOfftargetBtn.textContent = `Check severe off-targets (${rows.length})`;
    screenOfftargetBtn.disabled = Boolean(reason);
    if (offtargetHint) {
        if (reason) {
            offtargetHint.textContent = reason;
            offtargetHint.style.display = "";
        } else if (currentJobMode === "sequence") {
            // No blocker: keep the sequence-mode caveat permanently visible
            // below the button so the on-target limitation is never missed.
            offtargetHint.textContent = SEQUENCE_OFFTARGET_NOTE;
            offtargetHint.style.display = "";
        } else {
            offtargetHint.textContent = "";
            offtargetHint.style.display = "none";
        }
    }
    if (reason) {
        screenOfftargetBtn.title = reason;
    } else {
        // Short always-available explainer: the search is deliberately
        // limited to the most severe off-targets.
        screenOfftargetBtn.title =
            "Runs on the currently filtered guides. Only the most severe " +
            "off-targets are enumerated (CCTOP score at or below the " +
            "operator's threshold); guides without hits here can still " +
            "have milder off-targets.";
    }
}

// POST /offtarget for the filtered subset and poll the child job.
// The button stays blocked (spinner) until the screening job reaches a
// terminal state — ``screening`` is cleared in pollScreenJob, not here.
async function runScreening() {
    if (screening || !currentJobId) return;
    const rows = collectScreeningRows(currentFilteredGuideRows);
    const reason = screeningBlockReason(rows.length);
    if (reason) {
        showToast("warning", `Cannot start off-target search: ${reason}.`);
        return;
    }
    const body = { job_id: currentJobId, rows };
    if (currentJobMode === "sequence") {
        // Sequence mode: the on-target site cannot be excluded (the input
        // sequence's location in the reference genome is unknown). The
        // caveat is permanently shown below the search button, so starting
        // the search implies acceptance; the flag also keeps the server-side
        // guard happy, and the report repeats the caveat after the run.
        body.acknowledge_missing_on_target = true;
        if (otGenomeSelect && otGenomeSelect.value) {
            body.genome_key = otGenomeSelect.value;
        }
        sequenceScreenAcknowledged = true;
    } else {
        sequenceScreenAcknowledged = false;
    }
    screening = true;
    updateScreenButton();
    try {
        const res = await fetch("/offtarget", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(body),
        });
        const data = await res.json();
        if (!data.ok) {
            showToast("error", data.error || "Off-target search could not be started.");
            screening = false;
            updateScreenButton();
            return;
        }
        screenJobId = data.job_id;
        screenRevokeToken = data.revoke_token;
        if (offtargetStatus) {
            showAlert(offtargetStatus, "info",
                `Off-target search running for ${data.n_guides} guides\u2026`);
        }
        logLine(`Off-target screening queued (${data.n_guides} guides, job ${screenJobId.slice(0, 8)}\u2026).`);
        pollScreenJob(screenJobId);
    } catch (err) {
        showToast("error", `Off-target search failed to start: ${err}`);
        screening = false;
        updateScreenButton();
    }
}

function pollScreenJob(jobId) {
    const check = async () => {
        try {
            const res = await fetch("/status?job_id=" + encodeURIComponent(jobId));
            const data = await res.json();
            if (data.status === "complete" || data.status === "error" || data.status === "cancelled") {
                // The screening job is done either way — release the button.
                screening = false;
                if (data.status === "complete") {
                    await fetchScreenResults(jobId);
                } else if (data.status === "cancelled") {
                    if (offtargetStatus) showAlert(offtargetStatus, "warning", "Off-target search was cancelled.");
                    showToast("info", "Off-target search was cancelled.");
                } else {
                    if (offtargetStatus) showAlert(offtargetStatus, "error", "Off-target search failed. Check the Logging tab.");
                    showToast("error", "Off-target search failed. Check logs for details.");
                }
                updateScreenButton();
                return;
            }
            setTimeout(check, 1000);
        } catch (err) {
            setTimeout(check, 3000);
        }
    };
    setTimeout(check, 800);
}

async function fetchScreenResults(jobId) {
    const res = await fetch("/results?job_id=" + encodeURIComponent(jobId));
    const data = await res.json();
    if (!data.ok || !data.records || data.records.length === 0) {
        if (offtargetStatus) showAlert(offtargetStatus, "error", data.message || "No off-target results returned.");
        return;
    }
    screenRecords = data.records;
    // Build the verdict lookup used by the Guides-table verdict column, the
    // "hide failing guides" filter and the toggle title.
    otVerdictMap = new Map();
    let nFail = 0;
    for (const rec of screenRecords) {
        const v = offtargetVerdict(rec);
        if (!v) continue;
        otVerdictMap.set(JSON.stringify(screenRowKey(rec)), v);
        if (v.verdict === "fail") nFail++;
    }
    if (otHideFailedToggle) {
        otHideFailedToggle.disabled = false;
        otHideFailedToggle.title = nFail > 0
            ? `${nFail} of ${screenRecords.length} screened guides failed the off-target check.`
            : "No screened guide failed the off-target check.";
    }
    if (offtargetPlaceholder) offtargetPlaceholder.style.display = "none";
    if (offtargetStatus) {
        if (sequenceScreenAcknowledged) {
            // Sequence-mode screening ran with the on-target caveat
            // acknowledged; repeat it so the report is read correctly.
            showAlert(offtargetStatus, "warning",
                (data.message || "Off-target screening complete.") +
                " Sequence-mode caveat: the on-target site could not be " +
                "excluded — a perfect-match hit (worst score 0.0) may be " +
                "the on-target itself rather than a true off-target.");
        } else {
            showAlert(offtargetStatus, "success", data.message || "Off-target screening complete.");
        }
    }
    renderOfftargetTable();
    if (downloadOfftargetResults && screenJobId) {
        downloadOfftargetResults.href =
            "/download/results?job_id=" + encodeURIComponent(screenJobId);
        downloadOfftargetResults.style.display = "inline-flex";
    }
    // Stay on the Guides tab: the verdict column just appeared there. Add a
    // status line next to the guide counts and toast the outcome; the
    // report tab keeps its own banner for when the user visits it.
    renderResultsTable();
    const summary = `Off-target check: ${screenRecords.length} screened, ` +
        `${nFail} fail, ${screenRecords.length - nFail} pass.`;
    if (resultsGuideStatus) {
        let line = document.getElementById("otCheckStatus");
        if (!line) {
            line = document.createElement("div");
            line.id = "otCheckStatus";
            line.className = "ot-check-status";
            resultsGuideStatus.appendChild(line);
        }
        line.textContent = summary;
    }
    showToast("info", summary);
}

// Columns of the compact off-target report: the verdict plus the
// screening-relevant columns only — not a copy of the Guides table. The
// full screened records remain available via the TSV download.
const OFFTARGET_REPORT_COLUMNS = [
    "offtarget_verdict",
    "guide_rna_sequence",
    "base_editor",
    "worst_offtarget_score",
    "offtarget_count",
    "worst_offtarget_mismatch_pos",
    "worst_offtarget_chr",
    "worst_offtarget_pos",
    "worst_offtarget_strand",
    "worst_offtarget_sequence",
];

function renderOfftargetTable() {
    const wrap = offtargetPanel.querySelector(".results-table-wrap");
    if (!wrap) {
        const div = document.createElement("div");
        div.className = "results-table-wrap";
        offtargetPanel.appendChild(div);
    }
    const target = offtargetPanel.querySelector(".results-table-wrap");
    if (screenRecords.length === 0) {
        target.innerHTML = `<p style="color:var(--muted);padding:16px;">No off-target results returned.</p>`;
        return;
    }
    // Show only the report columns that the records actually carry (the
    // verdict column is always first and always shown).
    const present = new Set(Object.keys(screenRecords[0]));
    const columns = OFFTARGET_REPORT_COLUMNS.filter(
        c => c === "offtarget_verdict" || present.has(c));
    const displayColumns = columns.map(columnLabel);
    target.innerHTML = renderTable(
        screenRecords,
        columns,
        displayColumns,
        null,
        (row, column, value) => {
            if (column === "offtarget_verdict") return formatReportVerdictCell(row);
            const seqProps = formatSequencePropsCell(row, column, value);
            if (seqProps !== null) return seqProps;
            const otCell = formatOfftargetCell(row, column, value);
            if (otCell !== null) return otCell;
            return formatDnaSequenceCell(row, column, value);
        },
        resultsColumnTooltips()
    );
}

// Toggle the Summary sub-tab placeholder text.
function setSummaryPlaceholder(visible, text, tone) {
    const ph = document.getElementById("summaryPlaceholder");
    if (!ph) return;
    if (visible) {
        ph.textContent = text || "Summary will appear here after the pipeline completes.";
        // "Nothing found" states (no results / no guides match) render in
        // dark red so they stand out; pending/idle states stay muted grey.
        ph.style.color = tone === "danger" ? "#b91c1c" : "var(--muted)";
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
        setSummaryPlaceholder(true, "No guides match the current filters.", "danger");
        return;
    }

    // Targets without a guide RNA are intentionally excluded from the
    // per-editor summary; they are still visible in the Guides table.
    const realRows = filtered.filter(r => !isNoGuideRow(r));
    if (realRows.length === 0) {
        destroySummaryChart();
        resultsSummary.innerHTML = "";
        resultsSummary.style.display = "none";
        setSummaryPlaceholder(true, "No guides match the current filters.", "danger");
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

// "Hide failing guides": re-render the guides table with the off-target
// verdict filter applied/removed. Only meaningful once a screening
// snapshot exists for the current job (the checkbox is disabled before).
if (otHideFailedToggle) {
    otHideFailedToggle.addEventListener("change", () => {
        otHideFailed = otHideFailedToggle.checked;
        if (currentRecords.length > 0) renderResultsTable();
    });
}
