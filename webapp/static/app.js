const state = {
  bootstrap: null,
  mode: "library",
  selectedDrug: null,
  selectedScenario: "healthy_reference",
  selectedCommunity: "",
  customSessionId: null,
  customProfile: null,
  toastTimer: null,
  latestTrajectory: [],
};

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function formatNumber(value, digits = 3) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) {
    return "N/A";
  }
  return Number(value).toFixed(digits);
}

function formatPercent(value, digits = 1) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) {
    return "N/A";
  }
  return `${(Number(value) * 100).toFixed(digits)}%`;
}

function labelClass(label) {
  const key = String(label || "").toLowerCase().replace(/_/g, "-");
  return `label-${key || "unknown"}`;
}

function showToast(message, isError = false) {
  const toast = document.getElementById("toast");
  toast.textContent = message;
  toast.style.background = isError ? "rgba(127, 29, 29, 0.94)" : "rgba(31, 47, 42, 0.92)";
  toast.classList.add("visible");
  if (state.toastTimer) {
    window.clearTimeout(state.toastTimer);
  }
  state.toastTimer = window.setTimeout(() => toast.classList.remove("visible"), 2800);
}

async function fetchJson(url, options = {}) {
  const response = await fetch(url, options);
  const payload = await response.json();
  if (!response.ok) {
    throw new Error(payload.error || "Request failed");
  }
  return payload;
}

function createStatCard(label, value, caption = "") {
  return `
    <div class="stat-card">
      <span class="label">${escapeHtml(label)}</span>
      <span class="value">${value}</span>
      <span class="caption">${escapeHtml(caption)}</span>
    </div>
  `;
}

function createMetaChip(label, value) {
  if (!value) {
    return "";
  }
  return `<span class="meta-chip"><strong>${escapeHtml(label)}</strong>&nbsp;${escapeHtml(value)}</span>`;
}

function createStatusPill(label) {
  if (!label) {
    return `<span class="meta-chip">N/A</span>`;
  }
  return `<span class="status-pill ${labelClass(label)}">${escapeHtml(label)}</span>`;
}

function isCustomMode() {
  return state.mode === "custom" && Boolean(state.customSessionId);
}

function currentMicrobeRows() {
  if (isCustomMode()) {
    return state.bootstrap?.custom_microbes || state.bootstrap?.microbes || [];
  }
  return state.bootstrap?.microbes || [];
}

function currentPanelSize() {
  const summary = state.bootstrap?.summary || {};
  if (isCustomMode()) {
    return Number(summary.n_custom_microbes || currentMicrobeRows().length || 0);
  }
  return Number(summary.n_microbes || currentMicrobeRows().length || 0);
}

function currentCommunityTablePath() {
  const manualInput = document.getElementById("communityTablePathInput");
  const presetSelect = document.getElementById("communityTableSelect");
  const manualPath = manualInput ? manualInput.value.trim() : "";
  if (manualPath) {
    return manualPath;
  }
  if (presetSelect && presetSelect.value) {
    return presetSelect.value;
  }
  return "";
}

function renderMicrobePanelNote() {
  const note = document.getElementById("microbePanelNote");
  const hint = document.getElementById("microbeSelectionHint");
  if (!note || !hint) {
    return;
  }

  const count = currentPanelSize();
  if (isCustomMode()) {
    note.innerHTML = `<strong>新药预测默认直接展示 ${escapeHtml(count)} 个微生物整面板结果。</strong>`;
    hint.textContent = "页面展示的是 83 菌整面板分布、Top 菌列表和汇总指标，无需手动选择单个菌。";
    return;
  }

  note.innerHTML = `<strong>库内药物默认直接展示 ${escapeHtml(count)} 个微生物整面板结果。</strong>`;
  hint.textContent = "库内药物模式会展示当前药物在整面板上的整体分布、Top 菌列表和汇总指标。";
}

function renderModeBanner() {
  const banner = document.getElementById("predictionModeBanner");
  const communityPath = currentCommunityTablePath();
  const communitySuffix = communityPath ? " · 真实 cohort 初始化" : "";
  if (isCustomMode()) {
    const name = state.customProfile?.drug?.chemical_name || "Custom drug";
    const count = currentPanelSize();
    banner.textContent = `当前模式：新药 SMILES 预测 · ${name} · ${count} 菌扩展面板${communitySuffix}`;
    return;
  }
  if (!state.selectedDrug) {
    banner.textContent = "当前模式：摘要浏览 · 等待选择库内药物或输入新药 SMILES";
    return;
  }
  const count = currentPanelSize();
  banner.textContent = `当前模式：库内药物查询 · ${count} 菌整面板${communitySuffix}`;
}

function renderHeroStats(summary) {
  document.getElementById("heroStats").innerHTML = [
    createStatCard("Pairs", summary.n_pairs, "整合 Step 1 + Step 2 可用预测对"),
    createStatCard("Drugs", summary.n_drugs, "当前网页可直接查询的库内药物"),
    createStatCard("Microbes", summary.n_microbes, "当前覆盖的菌株/物种面板"),
    createStatCard("Applicable", summary.n_applicable_pairs, "落在当前 applicability 范围内的 pair"),
  ].join("");
}

function renderBootstrapSummary(summary) {
  const step1 = summary.step1_counts || {};
  const step2 = summary.step2_counts || {};
  document.getElementById("bootstrapSummary").innerHTML = [
    createStatCard("Step 1 Inhibit", step1.inhibit ?? 0, "药物-微生物抑制预测数"),
    createStatCard("Step 1 Promote", step1.promote ?? 0, "药物-微生物促进预测数"),
    createStatCard("Step 2 Metabolized", step2.metabolized ?? 0, "预测会被代谢的 pair"),
    createStatCard("Step 2 Not Metabolized", step2.not_metabolized ?? 0, "预测不会被代谢的 pair"),
  ].join("");
}

function renderDemoRanking(rows) {
  const container = document.getElementById("demoRanking");
  if (!rows || !rows.length) {
    container.innerHTML = `<div class="empty-state">当前没有可展示的 Step 3 demo ranking。</div>`;
    return;
  }
  container.innerHTML = rows
    .map(
      (row, index) => `
        <div class="ranking-row">
          <div class="ranking-rank">${index + 1}</div>
          <div class="ranking-meta">
            <strong>${escapeHtml(row.chemical_name || row.prestwick_id)}</strong>
            <span>${escapeHtml(row.prestwick_id || "")} · ${escapeHtml(row.scenario_name || "scenario")}</span>
          </div>
          <div class="ranking-score">${formatNumber(row.development_score, 2)}</div>
        </div>
      `
    )
    .join("");
}

function populateDrugSelect() {
  const select = document.getElementById("drugSelect");
  const filterText = document.getElementById("drugFilterInput").value.trim().toLowerCase();
  const drugs = (state.bootstrap?.drugs || []).filter((row) => {
    if (!filterText) {
      return true;
    }
    return (
      String(row.chemical_name || "").toLowerCase().includes(filterText) ||
      String(row.prestwick_id || "").toLowerCase().includes(filterText)
    );
  });
  const hasSelectedDrug = drugs.some((drug) => drug.prestwick_id === state.selectedDrug);
  if (!hasSelectedDrug) {
    state.selectedDrug = null;
  }

  select.innerHTML = [
    `<option value="">请选择库内药物后再加载预测结果</option>`,
    ...drugs.map(
      (drug) =>
        `<option value="${escapeHtml(drug.prestwick_id)}">${escapeHtml(drug.chemical_name)} (${escapeHtml(
          drug.prestwick_id
        )})</option>`
    ),
  ].join("");
  select.value = state.selectedDrug || "";
}

function populateScenarioSelect() {
  const select = document.getElementById("scenarioSelect");
  const scenarios = state.bootstrap?.scenarios || [];
  select.innerHTML = scenarios
    .map(
      (scenario) =>
        `<option value="${escapeHtml(scenario.scenario_name)}">${escapeHtml(scenario.scenario_name)}</option>`
    )
    .join("");
  if (!scenarios.some((scenario) => scenario.scenario_name === state.selectedScenario)) {
    state.selectedScenario = scenarios[0]?.scenario_name || "healthy_reference";
  }
  select.value = state.selectedScenario;
}

function populateCommunitySelect() {
  const select = document.getElementById("communityTableSelect");
  const communities = state.bootstrap?.cohort_communities || [];
  select.innerHTML = [
    `<option value="">不使用真实 cohort，继续用内置场景</option>`,
    ...communities.map(
      (community) =>
        `<option value="${escapeHtml(community.community_table_path)}">${escapeHtml(community.label)}</option>`
    ),
  ].join("");
  select.value = state.selectedCommunity || "";
}

function syncCommunityModeUi() {
  const scenarioSelect = document.getElementById("scenarioSelect");
  const compareButton = document.getElementById("compareScenariosButton");
  const hint = document.getElementById("communityTableHint");
  const activeCommunity = currentCommunityTablePath();
  const usingCommunity = Boolean(activeCommunity);
  scenarioSelect.disabled = usingCommunity;
  compareButton.textContent = usingCommunity ? "查看当前 cohort 摘要" : "比较全部内置场景";
  if (hint) {
    hint.textContent = usingCommunity
      ? "当前已启用真实 community_table。Step 3 会以该群落初始化，场景下拉会暂时停用。"
      : "可直接填写 community table 路径。手动路径优先级高于上面的预设下拉。";
  }
}

function renderSelectedDrugMeta(drug) {
  const container = document.getElementById("selectedDrugMeta");
  if (!drug) {
    container.innerHTML = "";
    return;
  }

  const modeLabel = isCustomMode() ? "New SMILES" : "Library";
  const canonicalSmiles = String(drug.canonical_smiles_rdkit || drug.smiles || "").trim();
  const shortSmiles = canonicalSmiles.length > 36 ? `${canonicalSmiles.slice(0, 36)}...` : canonicalSmiles;

  container.innerHTML = [
    createMetaChip("Mode", modeLabel),
    createMetaChip("Drug", drug.chemical_name),
    createMetaChip("ID", drug.prestwick_id),
    createMetaChip("Class", drug.therapeutic_class),
    createMetaChip("Effect", drug.therapeutic_effect),
    createMetaChip("Formula", drug.molecular_formula),
    createMetaChip("Scaffold", drug.murcko_scaffold),
    createMetaChip("SMILES", shortSmiles),
  ].join("");
}

function renderBarList(containerId, rows, config) {
  const container = document.getElementById(containerId);
  if (!rows || !rows.length) {
    container.innerHTML = `<div class="empty-state">暂无可展示数据。</div>`;
    return;
  }

  if (config.mode === "diverging") {
    const maxAbs = Math.max(...rows.map((row) => Math.abs(Number(row[config.valueKey] || 0))), 0.001);
    container.innerHTML = rows
      .map((row) => {
        const value = Number(row[config.valueKey] || 0);
        const width = `${(Math.abs(value) / maxAbs) * 50}%`;
        const left = value < 0 ? `calc(50% - ${width})` : "50%";
        const color = value < 0 ? "#dc2626" : "#16a34a";
        return `
          <div class="bar-row">
            <div class="bar-label">${escapeHtml(row[config.labelKey] || "Unknown")}</div>
            <div class="diverging-track">
              <div class="diverging-fill" style="left:${left}; width:${width}; background:${color};"></div>
            </div>
            <div class="bar-value">${formatNumber(value, 3)}</div>
          </div>
        `;
      })
      .join("");
    return;
  }

  const maxValue = Math.max(...rows.map((row) => Number(row[config.valueKey] || 0)), 0.001);
  container.innerHTML = rows
    .map((row) => {
      const value = Number(row[config.valueKey] || 0);
      const width = `${(value / maxValue) * 100}%`;
      return `
        <div class="bar-row">
          <div class="bar-label">${escapeHtml(row[config.labelKey] || "Unknown")}</div>
          <div class="bar-track">
            <div class="bar-fill" style="width:${width}; background:${config.color};"></div>
          </div>
          <div class="bar-value">${formatNumber(value, 3)}</div>
        </div>
      `;
    })
    .join("");
}

function renderDetailList(containerId, rows) {
  const container = document.getElementById(containerId);
  container.innerHTML = rows
    .map(
      (row) => `
        <div class="detail-row">
          <span class="detail-term">${escapeHtml(row.label)}</span>
          <span class="detail-value">${row.value}</span>
        </div>
      `
    )
    .join("");
}

function renderTableBody(bodyId, rows, mapper, emptyColspan = 4) {
  const tbody = document.getElementById(bodyId);
  if (!rows || !rows.length) {
    tbody.innerHTML = `<tr><td colspan="${emptyColspan}">暂无数据</td></tr>`;
    return;
  }
  tbody.innerHTML = rows.map(mapper).join("");
}

function resetStep1Display() {
  document.getElementById("step1PairStats").innerHTML = [
    createStatCard("Panel Scope", "N/A", "选择库内药物或输入新药后显示"),
    createStatCard("Inhibit Pairs", "N/A", "等待预测结果"),
    createStatCard("Mean Inhibit Prob", "N/A", "等待预测结果"),
    createStatCard("Mean Effect Score", "N/A", "等待预测结果"),
  ].join("");
  document.getElementById("step1EffectBars").innerHTML = `<div class="empty-state">请选择库内药物或输入新药 SMILES。</div>`;
  renderDetailList("step1PairDetails", [
    { label: "Status", value: "等待加载预测结果" },
    { label: "Panel Scope", value: "N/A" },
    { label: "Inhibit", value: "N/A" },
    { label: "Promote", value: "N/A" },
    { label: "No Effect", value: "N/A" },
    { label: "Strongest Effect", value: "N/A" },
  ]);
  renderTableBody("step1TableBody", [], () => "", 4);
}

function renderStep1(profile) {
  const aggregated = profile?.aggregated || {};
  const counts = aggregated.step1_counts || {};
  const panelRows = profile.panel_effect_microbes || profile.top_effect_microbes || [];
  const strongestRow = panelRows[0] || null;
  const panelSize = currentPanelSize();
  const pairStats = document.getElementById("step1PairStats");
  pairStats.innerHTML = [
    createStatCard("Panel Scope", `${panelSize}`, "当前药物已在整面板微生物上完成预测"),
    createStatCard(
      "Inhibit Pairs",
      counts.inhibit ?? 0,
      "预测为抑制的微生物 pair 数"
    ),
    createStatCard(
      "Mean Inhibit Prob",
      formatNumber(aggregated.mean_predicted_inhibit_probability, 3),
      "面板平均抑制概率"
    ),
    createStatCard(
      "Mean Effect Score",
      formatNumber(aggregated.mean_predicted_effect_score, 3),
      "面板平均连续效应"
    ),
  ].join("");

  renderBarList("step1EffectBars", panelRows.slice(0, 8), {
    mode: "diverging",
    labelKey: "microbe_label",
    valueKey: "predicted_effect_score",
  });

  renderDetailList("step1PairDetails", [
    { label: "Panel Scope", value: `${escapeHtml(panelSize)} microbes` },
    { label: "Inhibit", value: escapeHtml(String(counts.inhibit ?? 0)) },
    { label: "Promote", value: escapeHtml(String(counts.promote ?? 0)) },
    { label: "No Effect", value: escapeHtml(String(counts.no_effect ?? 0)) },
    { label: "AMR 修正", value: escapeHtml(String(aggregated.amr_corrected_pairs ?? 0)) },
    {
      label: "Strongest Effect",
      value: strongestRow ? escapeHtml(strongestRow.microbe_label || strongestRow.nt_code || "N/A") : "N/A",
    },
    {
      label: "Top Effect Score",
      value: strongestRow ? formatNumber(strongestRow.predicted_effect_score, 3) : "N/A",
    },
  ]);

  renderTableBody(
    "step1TableBody",
    panelRows,
    (row) => {
      const amrNote = row.amr_correction_applied
        ? `<br /><span class="muted">AMR 修正: ${escapeHtml(row.amr_expected_phenotype || "resistant")} prior</span>`
        : "";
      const rawLabelNote =
        row.amr_correction_applied && row.raw_predicted_effect_label
          ? `<br /><span class="muted">原始: ${escapeHtml(row.raw_predicted_effect_label)}</span>`
          : "";
      const rawProbNote =
        row.amr_correction_applied && row.raw_predicted_inhibit_probability !== null && row.raw_predicted_inhibit_probability !== undefined
          ? `<br /><span class="muted">原始 ${formatNumber(row.raw_predicted_inhibit_probability, 3)}</span>`
          : "";
      const rawScoreNote =
        row.amr_correction_applied && row.raw_predicted_effect_score !== null && row.raw_predicted_effect_score !== undefined
          ? `<br /><span class="muted">原始 ${formatNumber(row.raw_predicted_effect_score, 3)}</span>`
          : "";
      return `
        <tr>
          <td><strong>${escapeHtml(row.microbe_label || row.nt_code)}</strong><br /><span class="muted">${escapeHtml(
            row.genus || ""
          )} · ${escapeHtml(row.phylum || "")}</span>${amrNote}</td>
          <td>${createStatusPill(row.predicted_effect_label)}${rawLabelNote}</td>
          <td>${formatNumber(row.predicted_inhibit_probability, 3)}${rawProbNote}</td>
          <td>${formatNumber(row.predicted_effect_score, 3)}${rawScoreNote}</td>
        </tr>
      `;
    },
    4
  );
}

function resetStep2Display() {
  document.getElementById("step2PairStats").innerHTML = [
    createStatCard("Panel Scope", "N/A", "选择库内药物或输入新药后显示"),
    createStatCard("Metabolized Pairs", "N/A", "等待预测结果"),
    createStatCard("Mean Metabolized Prob", "N/A", "等待预测结果"),
    createStatCard("Applicability Rate", "N/A", "等待预测结果"),
  ].join("");
  document.getElementById("step2MetabolismBars").innerHTML = `<div class="empty-state">请选择库内药物或输入新药 SMILES。</div>`;
  renderDetailList("step2PairDetails", [
    { label: "Status", value: "等待加载预测结果" },
    { label: "Panel Scope", value: "N/A" },
    { label: "Metabolized", value: "N/A" },
    { label: "Not Metabolized", value: "N/A" },
    { label: "Applicability Rate", value: "N/A" },
    { label: "Top Metabolizer", value: "N/A" },
  ]);
  renderTableBody("step2TableBody", [], () => "", 4);
}

function renderStep2(profile) {
  const aggregated = profile?.aggregated || {};
  const counts = aggregated.step2_counts || {};
  const panelRows = profile.panel_metabolism_microbes || profile.top_metabolism_microbes || [];
  const strongestRow = panelRows[0] || null;
  const panelSize = currentPanelSize();
  document.getElementById("step2PairStats").innerHTML = [
    createStatCard("Panel Scope", `${panelSize}`, "当前药物已在整面板微生物上完成预测"),
    createStatCard(
      "Metabolized Pairs",
      counts.metabolized ?? 0,
      "预测会代谢该药物的微生物 pair 数"
    ),
    createStatCard(
      "Mean Metabolized Prob",
      formatNumber(aggregated.mean_predicted_metabolized_probability, 3),
      "面板平均代谢概率"
    ),
    createStatCard(
      "Applicability Rate",
      formatPercent(aggregated.applicability_rate, 1),
      "落在当前 applicability 范围内的比例"
    ),
  ].join("");

  renderBarList("step2MetabolismBars", panelRows.slice(0, 8), {
    mode: "standard",
    labelKey: "microbe_label",
    valueKey: "predicted_metabolized_probability",
    color: "linear-gradient(135deg, #d97706 0%, #f59e0b 100%)",
  });

  renderDetailList("step2PairDetails", [
    { label: "Panel Scope", value: `${escapeHtml(panelSize)} microbes` },
    { label: "Metabolized", value: escapeHtml(String(counts.metabolized ?? 0)) },
    { label: "Not Metabolized", value: escapeHtml(String(counts.not_metabolized ?? 0)) },
    { label: "Applicability Rate", value: formatPercent(aggregated.applicability_rate, 1) },
    { label: "机制投影覆盖", value: formatPercent(aggregated.mechanism_projection_rate, 1) },
    { label: "反应类已投影", value: escapeHtml(String(aggregated.reaction_projection_pairs ?? 0)) },
    { label: "基因证据已投影", value: escapeHtml(String(aggregated.gene_projection_pairs ?? 0)) },
    {
      label: "Top Metabolizer",
      value: strongestRow ? escapeHtml(strongestRow.microbe_label || strongestRow.nt_code || "N/A") : "N/A",
    },
    {
      label: "Top Probability",
      value: strongestRow ? formatNumber(strongestRow.predicted_metabolized_probability, 3) : "N/A",
    },
    {
      label: "Top 反应类",
      value: strongestRow ? escapeHtml(strongestRow.predicted_reaction_class || "N/A") : "N/A",
    },
  ]);

  renderTableBody(
    "step2TableBody",
    panelRows,
    (row) => {
      const reactionNote = row.predicted_reaction_class
        ? `<br /><span class="muted">反应类: ${escapeHtml(row.predicted_reaction_class)}${
            row.predicted_reaction_confidence !== null && row.predicted_reaction_confidence !== undefined
              ? ` · 置信 ${formatNumber(row.predicted_reaction_confidence, 2)}`
              : ""
          }</span>`
        : "";
      const productNote =
        row.predicted_candidate_product_count && row.predicted_candidate_product_ids
          ? `<br /><span class="muted">候选产物: ${escapeHtml(row.predicted_candidate_product_ids)}</span>`
          : "";
      const geneNote =
        row.predicted_evidence_gene_count && row.predicted_evidence_gene_ids
          ? `<br /><span class="muted">基因证据: ${escapeHtml(row.predicted_evidence_gene_ids)}</span>`
          : "";
      return `
        <tr>
          <td><strong>${escapeHtml(row.microbe_label || row.nt_code)}</strong><br /><span class="muted">${escapeHtml(
            row.genus || ""
          )} · ${escapeHtml(row.phylum || "")}</span>${reactionNote}</td>
          <td>${createStatusPill(row.predicted_metabolism_label)}</td>
          <td>${formatNumber(row.predicted_metabolized_probability, 3)}${geneNote}</td>
          <td>${formatNumber(row.predicted_parent_depletion_fraction, 3)}${productNote}</td>
        </tr>
      `;
    },
    4
  );
}

function drawTrajectoryChart(points) {
  const canvas = document.getElementById("step3TrajectoryChart");
  const context = canvas.getContext("2d");
  const parent = canvas.parentElement;
  const width = Math.max(320, parent.clientWidth - 8);
  const height = 320;
  const ratio = window.devicePixelRatio || 1;
  canvas.width = width * ratio;
  canvas.height = height * ratio;
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;
  context.scale(ratio, ratio);

  context.clearRect(0, 0, width, height);

  if (!points || !points.length) {
    context.fillStyle = "#5f6f69";
    context.font = "14px Bahnschrift, Segoe UI, sans-serif";
    context.fillText("暂无轨迹数据", 16, 28);
    return;
  }

  const margin = { top: 20, right: 20, bottom: 32, left: 44 };
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const maxX = Math.max(...points.map((row) => Number(row.timepoint || 0)), 1);
  const minX = Math.min(...points.map((row) => Number(row.timepoint || 0)), 0);
  const yMin = 0;
  const yMax = 100;

  function xScale(value) {
    const ratioX = (value - minX) / Math.max(maxX - minX, 1);
    return margin.left + ratioX * innerWidth;
  }

  function yScale(value) {
    const ratioY = (value - yMin) / (yMax - yMin);
    return margin.top + innerHeight - ratioY * innerHeight;
  }

  context.strokeStyle = "rgba(95, 111, 105, 0.18)";
  context.lineWidth = 1;
  for (let i = 0; i <= 4; i += 1) {
    const yValue = yMin + ((yMax - yMin) / 4) * i;
    const y = yScale(yValue);
    context.beginPath();
    context.moveTo(margin.left, y);
    context.lineTo(width - margin.right, y);
    context.stroke();
    context.fillStyle = "#5f6f69";
    context.font = "12px Bahnschrift, Segoe UI, sans-serif";
    context.fillText(String(Math.round(yValue)), 8, y + 4);
  }

  context.strokeStyle = "rgba(31, 47, 42, 0.34)";
  context.beginPath();
  context.moveTo(margin.left, margin.top);
  context.lineTo(margin.left, height - margin.bottom);
  context.lineTo(width - margin.right, height - margin.bottom);
  context.stroke();

  const series = [
    {
      label: "Health Index",
      color: "#0f766e",
      values: points.map((row) => Number(row.health_index || 0)),
    },
    {
      label: "Parent Retention ×100",
      color: "#d97706",
      values: points.map((row) => Number(row.parent_retention_ratio || 0) * 100),
    },
    {
      label: "Development Score",
      color: "#dc2626",
      values: points.map((row) => Number(row.development_score || 0)),
    },
  ];

  series.forEach((seriesItem) => {
    context.strokeStyle = seriesItem.color;
    context.lineWidth = 3;
    context.beginPath();
    seriesItem.values.forEach((value, index) => {
      const x = xScale(Number(points[index].timepoint || 0));
      const y = yScale(value);
      if (index === 0) {
        context.moveTo(x, y);
      } else {
        context.lineTo(x, y);
      }
    });
    context.stroke();

    seriesItem.values.forEach((value, index) => {
      const x = xScale(Number(points[index].timepoint || 0));
      const y = yScale(value);
      context.fillStyle = "#ffffff";
      context.beginPath();
      context.arc(x, y, 4, 0, Math.PI * 2);
      context.fill();
      context.fillStyle = seriesItem.color;
      context.beginPath();
      context.arc(x, y, 2.6, 0, Math.PI * 2);
      context.fill();
    });
  });

  const legendY = 16;
  let legendX = margin.left;
  series.forEach((seriesItem) => {
    context.fillStyle = seriesItem.color;
    context.fillRect(legendX, legendY - 9, 12, 12);
    context.fillStyle = "#1f2f2a";
    context.font = "12px Bahnschrift, Segoe UI, sans-serif";
    context.fillText(seriesItem.label, legendX + 18, legendY + 1);
    legendX += context.measureText(seriesItem.label).width + 42;
  });
}

function renderScenarioGrid(rows) {
  const container = document.getElementById("step3ScenarioGrid");
  if (!rows || !rows.length) {
    container.innerHTML = `<div class="empty-state">先运行 Step 3，再比较场景。</div>`;
    return;
  }
  container.innerHTML = rows
    .map(
      (row) => `
        <div class="scenario-card">
          <strong>${escapeHtml(row.scenario_name)}</strong>
          <div class="muted">${escapeHtml(row.scenario_description || row.community_source || "")}</div>
          <div class="scenario-metrics">
            <span>健康终值: ${formatNumber(row.final_health_index, 2)}</span>
            <span>TCG签名分: ${formatNumber(row.final_tcg_health_index, 2)}</span>
            <span>TCG覆盖率: ${formatNumber(row.final_tcg_mapped_fraction, 3)}</span>
            <span>母药保留: ${formatNumber(row.final_parent_retention_ratio, 3)}</span>
            <span>新版总分: ${formatNumber(row.development_score, 2)}</span>
            <span>旧版总分: ${formatNumber(row.development_score_legacy, 2)}</span>
          </div>
        </div>
      `
    )
    .join("");
}

function renderStep3Breakdown(summary) {
  renderDetailList("step3ScoreBreakdown", [
    { label: "收益分", value: formatNumber(summary.benefit_subscore_final, 2) },
    { label: "风险分", value: formatNumber(summary.risk_subscore_final, 2) },
    { label: "母药保留收益", value: formatNumber(summary.efficacy_proxy_final, 2) },
    { label: "群落保真收益", value: formatNumber(summary.community_preservation_final, 2) },
    { label: "TCG签名分", value: formatNumber(summary.final_tcg_health_index, 2) },
    { label: "TCG Guild 1 占比", value: formatNumber(summary.final_tcg_guild_1_fraction, 3) },
    { label: "TCG Guild 2 占比", value: formatNumber(summary.final_tcg_guild_2_fraction, 3) },
    { label: "TCG覆盖率", value: formatNumber(summary.final_tcg_mapped_fraction, 3) },
    { label: "菌群失衡惩罚", value: formatNumber(summary.dysbiosis_penalty_final, 2) },
    { label: "适用域惩罚", value: formatNumber(summary.uncertainty_penalty_final, 2) },
    { label: "代谢负担惩罚", value: formatNumber(summary.metabolite_burden_penalty_final, 2) },
    { label: "旧版启发式分", value: formatNumber(summary.development_score_legacy, 2) },
  ]);
}

function renderStep3(result) {
  const summary = result.summary || {};
  state.latestTrajectory = result.trajectory_metrics || [];
  document.getElementById("step3TrajectoryChart").dataset.hasData = "true";
  document.getElementById("step3Summary").innerHTML = [
    createStatCard(
      "场景",
      escapeHtml(summary.scenario_name || "N/A"),
      summary.community_source ? `${summary.scenario_description || ""} · ${summary.community_source}` : summary.scenario_description || ""
    ),
    createStatCard("最终健康", formatNumber(summary.final_health_index, 2), "健康指数终值"),
    createStatCard(
      "TCG签名",
      formatNumber(summary.final_tcg_health_index, 2),
      `A core microbiome signature proxy，覆盖率 ${formatNumber(summary.final_tcg_mapped_fraction, 3)}`
    ),
    createStatCard("母药保留", formatNumber(summary.final_parent_retention_ratio, 3), "累计给药后剩余母药比例"),
    createStatCard("综合评分", formatNumber(summary.development_score, 2), "综合排序分，结合收益项与惩罚项"),
  ].join("");
  renderStep3Breakdown(summary);

  drawTrajectoryChart(result.trajectory_metrics || []);
  renderTableBody(
    "step3TableBody",
    (result.top_microbe_changes || []).slice(0, 10),
    (row) => `
      <tr>
        <td><strong>${escapeHtml(row.species_label || row.nt_code)}</strong></td>
        <td>${formatNumber(row.initial_abundance, 4)}</td>
        <td>${formatNumber(row.final_abundance, 4)}</td>
        <td>${formatNumber(row.delta_abundance, 4)}</td>
        <td>${formatNumber(row.fold_change, 3)}</td>
      </tr>
    `,
    5
  );
}

function resetStep3Display() {
  document.getElementById("step3Summary").innerHTML = [
    createStatCard("场景", "N/A", "运行 Step 3 后显示"),
    createStatCard("最终健康", "N/A", "健康指数终值"),
    createStatCard("TCG签名", "N/A", "等待 guild 映射后显示"),
    createStatCard("母药保留", "N/A", "累计给药后剩余母药比例"),
    createStatCard("综合评分", "N/A", "综合排序分"),
  ].join("");
  renderDetailList("step3ScoreBreakdown", [
    { label: "收益分", value: "N/A" },
    { label: "风险分", value: "N/A" },
    { label: "母药保留收益", value: "N/A" },
    { label: "群落保真收益", value: "N/A" },
    { label: "TCG签名分", value: "N/A" },
    { label: "TCG Guild 1 占比", value: "N/A" },
    { label: "TCG Guild 2 占比", value: "N/A" },
    { label: "TCG覆盖率", value: "N/A" },
    { label: "菌群失衡惩罚", value: "N/A" },
    { label: "适用域惩罚", value: "N/A" },
    { label: "代谢负担惩罚", value: "N/A" },
    { label: "旧版启发式分", value: "N/A" },
  ]);
  state.latestTrajectory = [];
  document.getElementById("step3TrajectoryChart").dataset.hasData = "false";
  drawTrajectoryChart([]);
  renderScenarioGrid([]);
  renderTableBody("step3TableBody", [], () => "", 5);
}

function focusStep3Panel() {
  const panel = document.getElementById("step3Panel");
  if (!panel) {
    return;
  }
  panel.scrollIntoView({ behavior: "smooth", block: "start" });
}

function renderStep3Loading(message = "正在运行 Step 3...") {
  document.getElementById("step3Summary").innerHTML = [
    createStatCard("场景", "运行中", message),
    createStatCard("最终健康", "...", "请稍等"),
    createStatCard("TCG签名", "...", "正在匹配 guild proxy"),
    createStatCard("母药保留", "...", "请稍等"),
    createStatCard("综合评分", "...", "请稍等"),
  ].join("");
  renderDetailList("step3ScoreBreakdown", [
    { label: "收益分", value: "..." },
    { label: "风险分", value: "..." },
    { label: "母药保留收益", value: "..." },
    { label: "群落保真收益", value: "..." },
    { label: "TCG签名分", value: "..." },
    { label: "TCG Guild 1 占比", value: "..." },
    { label: "TCG Guild 2 占比", value: "..." },
    { label: "TCG覆盖率", value: "..." },
    { label: "菌群失衡惩罚", value: "..." },
    { label: "适用域惩罚", value: "..." },
    { label: "代谢负担惩罚", value: "..." },
    { label: "旧版启发式分", value: "..." },
  ]);
}

function buildSimulationPayload() {
  return {
    scenario: document.getElementById("scenarioSelect").value,
    community_table_path: currentCommunityTablePath() || null,
    n_steps: Number(document.getElementById("nStepsInput").value),
    initial_dose: Number(document.getElementById("initialDoseInput").value),
    repeat_dose: Number(document.getElementById("repeatDoseInput").value),
    dosing_interval: Number(document.getElementById("dosingIntervalInput").value),
    drug_clearance_rate: Number(document.getElementById("drugClearanceInput").value),
    product_clearance_rate: Number(document.getElementById("productClearanceInput").value),
    metabolism_scale: Number(document.getElementById("metabolismScaleInput").value),
    effect_scale: Number(document.getElementById("effectScaleInput").value),
    ecology_strength: Number(document.getElementById("ecologyStrengthInput").value),
  };
}

function setButtonBusy(buttonId, isBusy, busyText = "处理中...") {
  const button = document.getElementById(buttonId);
  if (!button.dataset.defaultText) {
    button.dataset.defaultText = button.textContent;
  }
  button.disabled = isBusy;
  button.style.opacity = isBusy ? "0.72" : "1";
  button.textContent = isBusy ? busyText : button.dataset.defaultText;
}

function clearCustomState() {
  state.mode = "library";
  state.customSessionId = null;
  state.customProfile = null;
}

function resetPredictionDisplays() {
  renderSelectedDrugMeta(null);
  resetStep1Display();
  resetStep2Display();
  resetStep3Display();
}

async function activateLibraryMode(showMessage = false) {
  clearCustomState();
  renderMicrobePanelNote();
  renderModeBanner();
  if (state.selectedDrug) {
    await loadPredictions();
    await runStep3Simulation();
  } else {
    resetPredictionDisplays();
  }
  if (showMessage) {
    showToast("已切回库内药物模式");
  }
}

async function loadPredictions() {
  try {
    if (isCustomMode()) {
      renderModeBanner();
      renderMicrobePanelNote();
      renderSelectedDrugMeta(state.customProfile?.drug);
      renderStep1(state.customProfile);
      renderStep2(state.customProfile);
      return;
    }

    if (!state.selectedDrug) {
      renderModeBanner();
      renderMicrobePanelNote();
      resetPredictionDisplays();
      return;
    }
    const profile = await fetchJson(`/api/drug-profile?drug=${encodeURIComponent(state.selectedDrug)}`);
    renderModeBanner();
    renderMicrobePanelNote();
    renderSelectedDrugMeta(profile.drug);
    renderStep1(profile);
    renderStep2(profile);
  } catch (error) {
    showToast(error.message, true);
  }
}

async function runStep3Simulation() {
  if (!isCustomMode() && !state.selectedDrug) {
    showToast("请先选择库内药物或输入新药 SMILES。", true);
    resetStep3Display();
    return false;
  }
  try {
    setButtonBusy("runSimulationButton", true, "模拟中...");
    focusStep3Panel();
    renderStep3Loading("正在计算单场景模拟结果");
    showToast(currentCommunityTablePath() ? "正在运行真实 cohort 初始化的 Step 3 模拟..." : "正在运行 Step 3 模拟...");
    const payload = buildSimulationPayload();
    let result;
    if (isCustomMode()) {
      result = await fetchJson("/api/custom-drug/step3/simulate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: state.customSessionId, ...payload }),
      });
    } else {
      if (!state.selectedDrug) {
        return;
      }
      result = await fetchJson("/api/step3/simulate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ drug: state.selectedDrug, ...payload }),
      });
    }
    renderStep3(result);
    showToast("Step 3 模拟完成");
    return true;
  } catch (error) {
    showToast(error.message, true);
    return false;
  } finally {
    setButtonBusy("runSimulationButton", false);
  }
}

async function loadScenarioGrid() {
  if (!isCustomMode() && !state.selectedDrug) {
    showToast("请先选择库内药物或输入新药 SMILES。", true);
    renderScenarioGrid([]);
    return;
  }
  try {
    setButtonBusy("compareScenariosButton", true, "对比中...");
    focusStep3Panel();
    const payload = buildSimulationPayload();
    let result;
    if (isCustomMode()) {
      result = await fetchJson("/api/custom-drug/step3/scenario-grid", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: state.customSessionId, ...payload }),
      });
    } else {
      if (!state.selectedDrug) {
        return;
      }
      result = await fetchJson("/api/step3/scenario-grid", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ drug: state.selectedDrug, ...payload }),
      });
    }
    renderScenarioGrid(result.scenario_summaries || []);
    showToast(currentCommunityTablePath() ? "已更新当前 cohort 摘要" : "已更新全部场景对比");
  } catch (error) {
    showToast(error.message, true);
  } finally {
    setButtonBusy("compareScenariosButton", false);
  }
}

async function predictCustomDrug() {
  const smiles = document.getElementById("customSmilesInput").value.trim();
  const drugName = document.getElementById("customDrugNameInput").value.trim();

  if (!smiles) {
    showToast("请输入新药 SMILES。", true);
    return;
  }

  try {
    setButtonBusy("predictCustomButton", true, "预测中...");
    showToast("正在根据 SMILES 生成 Step 1 / Step 2 预测...");
    const result = await fetchJson("/api/custom-drug/predict", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        smiles,
        drug_name: drugName || null,
      }),
    });

    state.mode = "custom";
    state.customSessionId = result.session_id;
    state.customProfile = result.profile;
    renderMicrobePanelNote();
    renderModeBanner();
    renderSelectedDrugMeta(result.profile.drug);
    renderStep1(result.profile);
    renderStep2(result.profile);
    resetStep3Display();
    showToast("新药 Step 1 / Step 2 预测完成");
    const step3Ok = await runStep3Simulation();
    if (step3Ok) {
      await loadScenarioGrid();
    }
  } catch (error) {
    showToast(error.message, true);
  } finally {
    setButtonBusy("predictCustomButton", false);
  }
}

function bindEvents() {
  document.getElementById("predictCustomButton").addEventListener("click", predictCustomDrug);
  document.getElementById("resetLibraryButton").addEventListener("click", () => activateLibraryMode(true));

  document.getElementById("customSmilesInput").addEventListener("keydown", (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
      predictCustomDrug();
    }
  });

  document.getElementById("drugFilterInput").addEventListener("input", async () => {
    const previous = state.selectedDrug;
    populateDrugSelect();
    if (state.mode === "library" && state.selectedDrug !== previous) {
      await loadPredictions();
      if (state.selectedDrug) {
        await runStep3Simulation();
      } else {
        resetStep3Display();
      }
    }
  });

  document.getElementById("drugSelect").addEventListener("change", async (event) => {
    state.selectedDrug = event.target.value;
    if (isCustomMode()) {
      clearCustomState();
      renderMicrobePanelNote();
      renderModeBanner();
      showToast("已切回库内药物模式");
    }
    await loadPredictions();
    if (state.selectedDrug) {
      await runStep3Simulation();
    } else {
      resetStep3Display();
    }
  });

  document.getElementById("scenarioSelect").addEventListener("change", async (event) => {
    state.selectedScenario = event.target.value;
    await runStep3Simulation();
  });

  document.getElementById("communityTableSelect").addEventListener("change", async (event) => {
    state.selectedCommunity = event.target.value;
    syncCommunityModeUi();
    renderModeBanner();
    if (isCustomMode() || state.selectedDrug) {
      await runStep3Simulation();
      renderScenarioGrid([]);
    }
  });

  document.getElementById("communityTablePathInput").addEventListener("change", async () => {
    syncCommunityModeUi();
    renderModeBanner();
    if (isCustomMode() || state.selectedDrug) {
      await runStep3Simulation();
      renderScenarioGrid([]);
    }
  });

  [
    "nStepsInput",
    "initialDoseInput",
    "repeatDoseInput",
    "dosingIntervalInput",
    "drugClearanceInput",
    "productClearanceInput",
    "metabolismScaleInput",
    "effectScaleInput",
    "ecologyStrengthInput",
  ].forEach((id) => {
    document.getElementById(id).addEventListener("change", async () => {
      const ok = await runStep3Simulation();
      if (ok) {
        renderScenarioGrid([]);
        showToast("Step 3 已按新参数重算；场景对比已清空，请按需重新比较。");
      }
    });
  });

  document.getElementById("refreshPredictionsButton").addEventListener("click", loadPredictions);
  document.getElementById("runSimulationButton").addEventListener("click", runStep3Simulation);
  document.getElementById("compareScenariosButton").addEventListener("click", loadScenarioGrid);

  window.addEventListener("resize", () => {
    const canvas = document.getElementById("step3TrajectoryChart");
    if (canvas.dataset.hasData === "true" && state.latestTrajectory) {
      drawTrajectoryChart(state.latestTrajectory);
    }
  });
}

async function bootstrap() {
  try {
    state.bootstrap = await fetchJson("/api/bootstrap");
    renderHeroStats(state.bootstrap.summary);
    renderBootstrapSummary(state.bootstrap.summary);
    renderDemoRanking(state.bootstrap.demo_candidates);

    state.selectedDrug = null;
    populateDrugSelect();
    renderMicrobePanelNote();
    populateScenarioSelect();
    populateCommunitySelect();
    syncCommunityModeUi();
    renderModeBanner();
    resetPredictionDisplays();
    bindEvents();
  } catch (error) {
    showToast(error.message, true);
  }
}

bootstrap();
