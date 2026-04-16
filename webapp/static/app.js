const state = {
  bootstrap: null,
  mode: "library",
  selectedDrug: null,
  selectedScenario: "healthy_reference",
  selectedDisease: "",
  selectedCommunity: "",
  customSessionId: null,
  customProfile: null,
  customSelectedPair: null,
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

function isMissingLike(value) {
  if (value === null || value === undefined) {
    return true;
  }
  if (typeof value === "number" && Number.isNaN(value)) {
    return true;
  }
  const text = String(value).trim().toLowerCase();
  return text === "" || text === "nan" || text === "none" || text === "null" || text === "n/a" || text === "na";
}

function textOrNA(value, fallback = "N/A") {
  return isMissingLike(value) ? fallback : String(value).trim();
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
  toast.setAttribute("aria-live", isError ? "assertive" : "polite");
  toast.classList.add("visible");
  if (state.toastTimer) {
    window.clearTimeout(state.toastTimer);
  }
  const duration = isError ? 9000 : 2800;
  state.toastTimer = window.setTimeout(() => toast.classList.remove("visible"), duration);
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

function createContextTag(label, tone = "neutral") {
  return `<span class="context-tag context-tag-${escapeHtml(tone)}">${escapeHtml(label)}</span>`;
}

function createEvidenceDetails(items) {
  const validItems = (items || []).filter((item) => String(item || "").trim());
  if (!validItems.length) {
    return "";
  }
  return `
    <details class="row-evidence">
      <summary>查看证据</summary>
      <div class="row-evidence-content">
        ${validItems.map((item) => `<div class="row-evidence-item">${escapeHtml(item)}</div>`).join("")}
      </div>
    </details>
  `;
}

function toFiniteNumber(value) {
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
}

function firstNonEmpty(rows, field, fallback = "") {
  for (const row of rows || []) {
    const raw = row?.[field];
    if (isMissingLike(raw)) {
      continue;
    }
    const value = String(raw).trim();
    if (value) {
      return value;
    }
  }
  return fallback;
}

function majorityLabel(rows, field, fallback = "") {
  const counts = new Map();
  for (const row of rows || []) {
    const raw = row?.[field];
    if (isMissingLike(raw)) {
      continue;
    }
    const value = String(raw).trim();
    if (!value) {
      continue;
    }
    counts.set(value, (counts.get(value) || 0) + 1);
  }
  if (!counts.size) {
    return fallback;
  }
  let best = fallback;
  let maxCount = -1;
  for (const [label, count] of counts.entries()) {
    if (count > maxCount) {
      best = label;
      maxCount = count;
    }
  }
  return best;
}

function meanField(rows, field) {
  const values = (rows || [])
    .map((row) => toFiniteNumber(row?.[field]))
    .filter((value) => value !== null);
  if (!values.length) {
    return null;
  }
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function anyTrue(rows, field) {
  return (rows || []).some((row) => Boolean(row?.[field]));
}

function uniqueJoinedValues(rows, field, maxItems = 3) {
  const values = [];
  const seen = new Set();
  for (const row of rows || []) {
    const rawValue = row?.[field];
    if (isMissingLike(rawValue)) {
      continue;
    }
    const raw = String(rawValue).trim();
    if (!raw || seen.has(raw)) {
      continue;
    }
    seen.add(raw);
    values.push(raw);
    if (values.length >= maxItems) {
      break;
    }
  }
  return values.join(" / ");
}

function aggregateMicrobeRows(rows, mode = "step1") {
  if (!rows || !rows.length) {
    return [];
  }
  const groups = new Map();
  for (const row of rows) {
    const displayName = String(row?.microbe_label || row?.species_label || row?.nt_code || "Unknown").trim();
    const key = displayName.toLowerCase();
    if (!groups.has(key)) {
      groups.set(key, { displayName, rows: [] });
    }
    groups.get(key).rows.push(row);
  }

  const aggregatedRows = Array.from(groups.values()).map((group) => {
    const groupRows = group.rows;
    const first = groupRows[0] || {};
    const ntCodes = [...new Set(groupRows.map((row) => String(row?.nt_code || "").trim()).filter(Boolean))];
    const base = {
      ...first,
      microbe_label: group.displayName,
      strain_count: ntCodes.length,
      strain_nt_codes: ntCodes,
      genus: firstNonEmpty(groupRows, "genus"),
      phylum: firstNonEmpty(groupRows, "phylum"),
      family: firstNonEmpty(groupRows, "family"),
      species_label: firstNonEmpty(groupRows, "species_label", group.displayName),
    };

    if (mode === "step2") {
      return {
        ...base,
        predicted_metabolism_label: majorityLabel(groupRows, "predicted_metabolism_label", first.predicted_metabolism_label),
        predicted_metabolized_probability: meanField(groupRows, "predicted_metabolized_probability"),
        predicted_parent_depletion_fraction: meanField(groupRows, "predicted_parent_depletion_fraction"),
        predicted_reaction_class: uniqueJoinedValues(groupRows, "predicted_reaction_class", 2),
        predicted_reaction_confidence: meanField(groupRows, "predicted_reaction_confidence"),
        predicted_enzyme_prior_flag: anyTrue(groupRows, "predicted_enzyme_prior_flag"),
        predicted_enzyme_support_score: meanField(groupRows, "predicted_enzyme_support_score"),
        predicted_enzyme_names: uniqueJoinedValues(groupRows, "predicted_enzyme_names", 2),
        predicted_enzyme_ids: uniqueJoinedValues(groupRows, "predicted_enzyme_ids", 2),
        predicted_enzyme_reaction_classes: uniqueJoinedValues(groupRows, "predicted_enzyme_reaction_classes", 2),
        predicted_enzyme_bond_targets: uniqueJoinedValues(groupRows, "predicted_enzyme_bond_targets", 2),
        predicted_evidence_gene_ids: uniqueJoinedValues(groupRows, "predicted_evidence_gene_ids", 2),
        predicted_candidate_product_ids: uniqueJoinedValues(groupRows, "predicted_candidate_product_ids", 2),
        predicted_evidence_gene_count: (groupRows || [])
          .map((row) => toFiniteNumber(row?.predicted_evidence_gene_count) || 0)
          .reduce((sum, value) => sum + value, 0),
        predicted_candidate_product_count: (groupRows || [])
          .map((row) => toFiniteNumber(row?.predicted_candidate_product_count) || 0)
          .reduce((sum, value) => sum + value, 0),
      };
    }

    return {
      ...base,
      predicted_effect_label: majorityLabel(groupRows, "predicted_effect_label", first.predicted_effect_label),
      predicted_inhibit_probability: meanField(groupRows, "predicted_inhibit_probability"),
      predicted_promote_probability_refined: meanField(groupRows, "predicted_promote_probability_refined"),
      predicted_effect_score: meanField(groupRows, "predicted_effect_score"),
      amr_correction_applied: anyTrue(groupRows, "amr_correction_applied"),
      amr_expected_phenotype: uniqueJoinedValues(groupRows, "amr_expected_phenotype", 1),
      raw_predicted_effect_label: majorityLabel(groupRows, "raw_predicted_effect_label", first.raw_predicted_effect_label),
      raw_predicted_inhibit_probability: meanField(groupRows, "raw_predicted_inhibit_probability"),
      raw_predicted_effect_score: meanField(groupRows, "raw_predicted_effect_score"),
      predicted_promote_support_type: majorityLabel(groupRows, "predicted_promote_support_type", first.predicted_promote_support_type),
      predicted_promote_evidence_type: majorityLabel(groupRows, "predicted_promote_evidence_type", first.predicted_promote_evidence_type),
      predicted_promote_support_score: meanField(groupRows, "predicted_promote_support_score"),
      predicted_cross_feeding_reference_flag: anyTrue(groupRows, "predicted_cross_feeding_reference_flag"),
      predicted_cross_feeding_support_microbe: uniqueJoinedValues(groupRows, "predicted_cross_feeding_support_microbe", 2),
      predicted_cross_feeding_match_mode: majorityLabel(groupRows, "predicted_cross_feeding_match_mode", first.predicted_cross_feeding_match_mode),
      predicted_cross_feeding_matched_term: uniqueJoinedValues(groupRows, "predicted_cross_feeding_matched_term", 2),
    };
  });

  if (mode === "step2") {
    return aggregatedRows.sort(
      (left, right) =>
        Math.abs(toFiniteNumber(right.predicted_metabolized_probability) || 0) -
        Math.abs(toFiniteNumber(left.predicted_metabolized_probability) || 0)
    );
  }

  return aggregatedRows.sort(
    (left, right) =>
      Math.abs(toFiniteNumber(right.predicted_effect_score) || 0) -
      Math.abs(toFiniteNumber(left.predicted_effect_score) || 0)
  );
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

function currentDiseaseName() {
  const select = document.getElementById("diseaseSelect");
  return select ? String(select.value || "").trim() : "";
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
  const diseaseName = currentDiseaseName();
  const communitySuffix = communityPath ? " · 真实 cohort 初始化" : "";
  const diseaseSuffix = diseaseName ? ` · 疾病背景 ${diseaseName}` : "";
  if (isCustomMode()) {
    const name = state.customProfile?.drug?.chemical_name || "Custom drug";
    const count = currentPanelSize();
    banner.textContent = `当前模式：新药 SMILES 预测 · ${name} · ${count} 菌扩展面板${diseaseSuffix}${communitySuffix}`;
    return;
  }
  if (!state.selectedDrug) {
    banner.textContent = "当前模式：摘要浏览 · 等待选择库内药物或输入新药 SMILES";
    return;
  }
  const count = currentPanelSize();
  banner.textContent = `当前模式：库内药物查询 · ${count} 菌整面板${diseaseSuffix}${communitySuffix}`;
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

function populateDiseaseSelect() {
  const select = document.getElementById("diseaseSelect");
  const diseases = state.bootstrap?.diseases || [];
  select.innerHTML = [
    `<option value="">不指定疾病背景，继续使用默认群落</option>`,
    ...diseases.map(
      (disease) =>
        `<option value="${escapeHtml(disease.disease_name)}">${escapeHtml(disease.disease_name)} (${escapeHtml(
          String(disease.microbe_relation_count || 0)
        )} 条菌群关系)</option>`
    ),
  ].join("");
  select.value = state.selectedDisease || "";
}

function syncCommunityModeUi() {
  const scenarioSelect = document.getElementById("scenarioSelect");
  const compareButton = document.getElementById("compareScenariosButton");
  const hint = document.getElementById("communityTableHint");
  const diseaseHint = document.getElementById("diseaseHint");
  const activeCommunity = currentCommunityTablePath();
  const usingCommunity = Boolean(activeCommunity);
  scenarioSelect.disabled = usingCommunity;
  compareButton.textContent = usingCommunity ? "查看当前 cohort 摘要" : "比较全部内置场景";
  if (hint) {
    hint.textContent = usingCommunity
      ? "当前已启用真实 community_table。Step 3 会以该群落初始化，场景下拉会暂时停用。"
      : "可直接填写 community table 路径。手动路径优先级高于上面的预设下拉。";
  }
  if (diseaseHint) {
    diseaseHint.textContent = usingCommunity
      ? "当前真实 community_table 优先，疾病背景只保留在解释层，不再覆盖 Step 3 初始化。"
      : "选择后会按疾病知识字典生成 Step 3 初始群落，并保留候选疾病解释。";
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

function renderDiseaseExplainCard(profile) {
  const selectedDisease = currentDiseaseName();
  const candidates = profile?.candidate_diseases || [];
  const marketed = profile?.marketed_disease_context || [];

  if (!profile) {
    renderDetailList("diseaseExplainCard", [
      { label: "状态", value: "等待加载药物结果" },
      { label: "说明", value: "当前还没有可展示的疾病解释。" },
    ]);
    return;
  }

  const selectedMatch =
    candidates.find((item) => String(item.disease_name || "").trim() === selectedDisease) ||
    candidates[0] ||
    null;
  const marketedLead = selectedDisease
    ? marketed.find((item) => String(item.disease_name || "").trim() === selectedDisease) || marketed[0] || null
    : marketed[0] || null;
  const evidence = selectedMatch?.evidence_examples?.[0] || null;

  renderDetailList("diseaseExplainCard", [
    { label: "当前疾病背景", value: selectedDisease ? escapeHtml(selectedDisease) : "未指定" },
    {
      label: "候选疾病",
      value: selectedMatch ? escapeHtml(selectedMatch.disease_name || "N/A") : "暂无候选疾病",
    },
    {
      label: "支持分",
      value: selectedMatch ? formatNumber(selectedMatch.support_score, 3) : "N/A",
    },
    {
      label: "匹配关系数",
      value: selectedMatch ? escapeHtml(String(selectedMatch.matched_relation_count || 0)) : "N/A",
    },
    {
      label: "代表证据",
      value: evidence
        ? `${escapeHtml(evidence.microbe_name_raw || evidence.matched_microbe || "N/A")} · 期望 ${escapeHtml(
            evidence.desired_step1_effect || "N/A"
          )} · 当前 ${escapeHtml(evidence.matched_effect_label || "N/A")}`
        : "暂无代表证据",
    },
    {
      label: "上市药物参考",
      value: marketedLead?.matched_market_drugs?.length
        ? escapeHtml(marketedLead.matched_market_drugs.join(" / "))
        : "暂无直接药物命中",
    },
    {
      label: "说明",
      value: selectedDisease
        ? "若未提供真实 community_table，Step 3 会优先用该疾病背景生成起始群落。"
        : "这里展示的是根据菌群方向性推断出的潜在疾病关联，不等同于临床适应证。",
    },
  ]);
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
    { label: "Mean Refined Promote Prob", value: "N/A" },
    { label: "Metabolism-supported Promote", value: "N/A" },
    { label: "No Effect", value: "N/A" },
    { label: "Strongest Effect", value: "N/A" },
  ]);
  renderTableBody("step1TableBody", [], () => "", 5);
  renderMechanismExplainCard(null, null);
  renderDiseaseExplainCard(null);
}

function describePromoteSupportType(value) {
  const mapping = {
    self_metabolism_supported: "自身代谢支持",
    self_metabolism_consistent: "自身代谢一致",
    cross_feeding_reference: "交叉喂养文献支持",
    weak_or_unspecified: "直接/宿主样或证据较弱",
  };
  return mapping[value] || value || "";
}

function describePromoteEvidenceType(value) {
  const mapping = {
    self_metabolism_supported_promote: "自身代谢促进",
    self_metabolism_consistent_promote: "自身代谢一致促进",
    cross_feeding_supported_promote: "交叉喂养促进",
    direct_or_host_like_promote: "直接或宿主环境促进",
    unspecified_promote: "未明确机制",
  };
  return mapping[value] || value || "";
}

function describeCrossFeedingMatchMode(value) {
  const mapping = {
    exact_compound: "精确化合物",
    compound_alias: "别名命中",
    keyword_family: "关键词家族",
    compound_family: "家族命中",
  };
  return mapping[value] || value || "";
}

function pickMechanismLead(profile) {
  const effectRows = profile?.panel_effect_microbes || profile?.top_effect_microbes || [];
  const metabolismRows = profile?.panel_metabolism_microbes || profile?.top_metabolism_microbes || [];
  if (!effectRows.length) {
    return null;
  }
  const metabolismMap = new Map(metabolismRows.map((row) => [row.nt_code, row]));
  const annotated = effectRows.map((row) => ({
    ...row,
    metabolism: metabolismMap.get(row.nt_code) || null,
  }));
  return (
    annotated.find((row) => Boolean(row.metabolism?.predicted_enzyme_prior_flag)) ||
    annotated.find((row) => row.predicted_cross_feeding_reference_flag) ||
    annotated.find(
      (row) =>
        row.predicted_promote_support_type === "self_metabolism_supported" ||
        row.predicted_promote_support_type === "self_metabolism_consistent"
    ) ||
    annotated.find((row) => row.predicted_effect_label === "promote") ||
    annotated[0]
  );
}

function renderMechanismExplainCard(profile, selectedPair = null) {
  const lead =
    selectedPair
      ? {
          nt_code: selectedPair.microbe?.nt_code,
          microbe_label: selectedPair.microbe?.microbe_label,
          predicted_effect_label: selectedPair.step1?.predicted_effect_label,
          predicted_promote_probability_refined: selectedPair.step1?.predicted_promote_probability_refined,
          predicted_promote_support_type: selectedPair.step1?.predicted_promote_support_type,
          predicted_promote_evidence_type: selectedPair.step1?.predicted_promote_evidence_type,
          predicted_cross_feeding_reference_flag: selectedPair.step1?.predicted_cross_feeding_reference_flag,
          predicted_cross_feeding_support_microbe: selectedPair.step1?.predicted_cross_feeding_support_microbe,
          predicted_cross_feeding_match_mode: selectedPair.step1?.predicted_cross_feeding_match_mode,
          predicted_cross_feeding_matched_term: selectedPair.step1?.predicted_cross_feeding_matched_term,
          metabolism: selectedPair.step2 || null,
        }
      : pickMechanismLead(profile);

  if (!lead) {
    renderDetailList("mechanismExplainCard", [
      { label: "状态", value: "等待加载预测结果" },
      { label: "说明", value: "当前还没有可展示的机制解释。" },
    ]);
    return;
  }

  const metabolism = lead.metabolism || {};
  const reactionClass = textOrNA(metabolism.predicted_reaction_class, "");
  const enzymeReactionClass = textOrNA(metabolism.predicted_enzyme_reaction_classes, "");
  const reactionDisplay = reactionClass || enzymeReactionClass || "N/A";
  const profileKey = String(
    profile?.confidence_breakdown?.drug_profile || profile?.aggregated?.confidence_breakdown?.drug_profile || ""
  )
    .trim()
    .toLowerCase();
  const summaryParts = [];
  if (lead.predicted_cross_feeding_reference_flag) {
    summaryParts.push("命中交叉喂养参考");
  }
  if (lead.predicted_promote_support_type) {
    summaryParts.push(describePromoteSupportType(lead.predicted_promote_support_type));
  }
  if (lead.predicted_promote_evidence_type) {
    summaryParts.push(describePromoteEvidenceType(lead.predicted_promote_evidence_type));
  }
  if (metabolism.predicted_enzyme_prior_flag) {
    summaryParts.push("酶先验支持");
  }
  if (profileKey === "sulfonamide_antifolate" && !metabolism.predicted_enzyme_prior_flag) {
    summaryParts.push("该药物以靶点抑菌机制为主，非代谢酶主导");
  }

  renderDetailList("mechanismExplainCard", [
    { label: "重点菌", value: escapeHtml(textOrNA(lead.microbe_label || lead.nt_code, "N/A")) },
    { label: "当前标签", value: createStatusPill(lead.predicted_effect_label || "N/A") },
    { label: "解释摘要", value: summaryParts.length ? escapeHtml(summaryParts.join(" / ")) : "暂无明确机制线索" },
    { label: "Promote 概率", value: formatNumber(lead.predicted_promote_probability_refined, 3) },
    { label: "供体菌", value: escapeHtml(textOrNA(lead.predicted_cross_feeding_support_microbe, "N/A")) },
    { label: "命中方式", value: escapeHtml(describeCrossFeedingMatchMode(lead.predicted_cross_feeding_match_mode) || "N/A") },
    { label: "命中词", value: escapeHtml(textOrNA(lead.predicted_cross_feeding_matched_term, "N/A")) },
    { label: "代谢概率", value: formatNumber(metabolism.predicted_metabolized_probability, 3) },
    { label: "母药消耗", value: formatNumber(metabolism.predicted_parent_depletion_fraction, 3) },
    { label: "反应类", value: escapeHtml(reactionDisplay) },
    { label: "酶先验支持", value: metabolism.predicted_enzyme_prior_flag ? "是" : "否" },
    { label: "酶支持分", value: formatNumber(metabolism.predicted_enzyme_support_score, 3) },
    { label: "候选酶", value: escapeHtml(textOrNA(metabolism.predicted_enzyme_names || metabolism.predicted_enzyme_ids, "N/A")) },
    { label: "酶反应类", value: escapeHtml(textOrNA(metabolism.predicted_enzyme_reaction_classes, "N/A")) },
    { label: "键靶点", value: escapeHtml(textOrNA(metabolism.predicted_enzyme_bond_targets, "N/A")) },
  ]);
}

function renderStep1(profile) {
  const aggregated = profile?.aggregated || {};
  const counts = aggregated.step1_counts || {};
  const rawPanelRows = profile.panel_effect_microbes || profile.top_effect_microbes || [];
  const panelRows = aggregateMicrobeRows(rawPanelRows, "step1");
  const strongestRow = panelRows[0] || null;
  const panelSize = currentPanelSize();
  const displaySize = panelRows.length;
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
    { label: "Display Scope", value: `${escapeHtml(displaySize)} unique names` },
    { label: "Inhibit", value: escapeHtml(String(counts.inhibit ?? 0)) },
    { label: "Promote", value: escapeHtml(String(counts.promote ?? 0)) },
    {
      label: "Mean Refined Promote Prob",
      value: formatNumber(aggregated.mean_predicted_promote_probability_refined, 3),
    },
    {
      label: "Metabolism-supported Promote",
      value: escapeHtml(String(aggregated.metabolism_supported_promote_pairs ?? 0)),
    },
    {
      label: "Cross-feeding Supported",
      value: escapeHtml(String(aggregated.cross_feeding_supported_promote_pairs ?? 0)),
    },
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
      const taxonomy = [row.genus, row.phylum].filter((value) => String(value || "").trim()).join(" · ");
      const tags = [];
      const evidence = [];
      if (Number(row.strain_count || 1) > 1) {
        tags.push(createContextTag(`${row.strain_count} strains`, "neutral"));
        evidence.push(`NT codes: ${(row.strain_nt_codes || []).join(", ")}`);
      }

      if (row.amr_correction_applied) {
        tags.push(createContextTag("AMR 修正", "warn"));
        evidence.push(`AMR prior: ${row.amr_expected_phenotype || "resistant"}`);
        if (row.raw_predicted_effect_label) {
          evidence.push(`原始标签: ${row.raw_predicted_effect_label}`);
        }
        if (row.raw_predicted_inhibit_probability !== null && row.raw_predicted_inhibit_probability !== undefined) {
          evidence.push(`原始 Inhibit 概率: ${formatNumber(row.raw_predicted_inhibit_probability, 3)}`);
        }
        if (row.raw_predicted_effect_score !== null && row.raw_predicted_effect_score !== undefined) {
          evidence.push(`原始 Effect Score: ${formatNumber(row.raw_predicted_effect_score, 3)}`);
        }
      }

      if (row.predicted_effect_label === "promote" && row.predicted_promote_support_type) {
        tags.push(createContextTag("Promote 支持", "good"));
        evidence.push(
          `支持来源: ${describePromoteSupportType(row.predicted_promote_support_type)} / ${describePromoteEvidenceType(
            row.predicted_promote_evidence_type
          )}`
        );
        evidence.push(`支持分: ${formatNumber(row.predicted_promote_support_score, 2)}`);
        if (row.predicted_cross_feeding_support_microbe) {
          evidence.push(`供体菌: ${row.predicted_cross_feeding_support_microbe}`);
        }
        if (row.predicted_cross_feeding_match_mode) {
          evidence.push(`命中方式: ${describeCrossFeedingMatchMode(row.predicted_cross_feeding_match_mode)}`);
        }
        if (row.predicted_cross_feeding_matched_term) {
          evidence.push(`命中词: ${row.predicted_cross_feeding_matched_term}`);
        }
      }

      return `
        <tr>
          <td>
            <div class="microbe-cell">
              <strong>${escapeHtml(row.microbe_label || row.nt_code)}</strong>
              <span class="muted">${escapeHtml(taxonomy || "未标注分类")}</span>
              ${tags.length ? `<div class="context-tags">${tags.join("")}</div>` : ""}
              ${createEvidenceDetails(evidence)}
            </div>
          </td>
          <td>${createStatusPill(row.predicted_effect_label)}</td>
          <td>${formatNumber(row.predicted_inhibit_probability, 3)}</td>
          <td>${formatNumber(row.predicted_promote_probability_refined, 3)}</td>
          <td>${formatNumber(row.predicted_effect_score, 3)}</td>
        </tr>
      `;
    },
    5
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
  const rawPanelRows = profile.panel_metabolism_microbes || profile.top_metabolism_microbes || [];
  const panelRows = aggregateMicrobeRows(rawPanelRows, "step2");
  const strongestRow = panelRows[0] || null;
  const panelSize = currentPanelSize();
  const displaySize = panelRows.length;
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
    { label: "Display Scope", value: `${escapeHtml(displaySize)} unique names` },
    { label: "Metabolized", value: escapeHtml(String(counts.metabolized ?? 0)) },
    { label: "Not Metabolized", value: escapeHtml(String(counts.not_metabolized ?? 0)) },
    { label: "Applicability Rate", value: formatPercent(aggregated.applicability_rate, 1) },
    { label: "机制投影覆盖", value: formatPercent(aggregated.mechanism_projection_rate, 1) },
    { label: "酶先验覆盖", value: formatPercent(aggregated.enzyme_prior_support_rate, 1) },
    { label: "反应类已投影", value: escapeHtml(String(aggregated.reaction_projection_pairs ?? 0)) },
    { label: "基因证据已投影", value: escapeHtml(String(aggregated.gene_projection_pairs ?? 0)) },
    { label: "酶先验支持对数", value: escapeHtml(String(aggregated.enzyme_prior_supported_pairs ?? 0)) },
    { label: "平均酶支持分", value: formatNumber(aggregated.mean_enzyme_support_score, 3) },
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
      const taxonomy = [row.genus, row.phylum].filter((value) => String(value || "").trim()).join(" · ");
      const tags = [];
      const evidence = [];
      if (Number(row.strain_count || 1) > 1) {
        tags.push(createContextTag(`${row.strain_count} strains`, "neutral"));
        evidence.push(`NT codes: ${(row.strain_nt_codes || []).join(", ")}`);
      }

      if (row.predicted_reaction_class) {
        tags.push(createContextTag("反应类", "neutral"));
        evidence.push(
          `反应类: ${row.predicted_reaction_class}${
            row.predicted_reaction_confidence !== null && row.predicted_reaction_confidence !== undefined
              ? ` (置信 ${formatNumber(row.predicted_reaction_confidence, 2)})`
              : ""
          }`
        );
      }

      if (row.predicted_enzyme_prior_flag) {
        tags.push(createContextTag("酶先验", "good"));
        evidence.push(
          `酶先验: ${row.predicted_enzyme_names || row.predicted_enzyme_ids || "supported"} · 支持 ${formatNumber(
            row.predicted_enzyme_support_score,
            2
          )}`
        );
        if (row.predicted_enzyme_reaction_classes) {
          evidence.push(`酶反应类: ${row.predicted_enzyme_reaction_classes}`);
        }
        if (row.predicted_enzyme_bond_targets) {
          evidence.push(`键靶点: ${row.predicted_enzyme_bond_targets}`);
        }
      }

      if (row.predicted_evidence_gene_count && row.predicted_evidence_gene_ids) {
        tags.push(createContextTag("基因证据", "neutral"));
        evidence.push(`基因证据: ${row.predicted_evidence_gene_ids}`);
      }

      if (row.predicted_candidate_product_count && row.predicted_candidate_product_ids) {
        tags.push(createContextTag("候选产物", "neutral"));
        evidence.push(`候选产物: ${row.predicted_candidate_product_ids}`);
      }

      return `
        <tr>
          <td>
            <div class="microbe-cell">
              <strong>${escapeHtml(row.microbe_label || row.nt_code)}</strong>
              <span class="muted">${escapeHtml(taxonomy || "未标注分类")}</span>
              ${tags.length ? `<div class="context-tags">${tags.join("")}</div>` : ""}
              ${createEvidenceDetails(evidence)}
            </div>
          </td>
          <td>${createStatusPill(row.predicted_metabolism_label)}</td>
          <td>${formatNumber(row.predicted_metabolized_probability, 3)}</td>
          <td>${formatNumber(row.predicted_parent_depletion_fraction, 3)}</td>
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
    disease_name: currentDiseaseName() || null,
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
  if (!button) {
    return;
  }
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
  state.customSelectedPair = null;
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
      renderMechanismExplainCard(state.customProfile, state.customSelectedPair);
      renderDiseaseExplainCard(state.customProfile);
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
    renderMechanismExplainCard(profile, null);
    renderDiseaseExplainCard(profile);
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
    state.customSelectedPair = result.selected_pair;
    renderMicrobePanelNote();
    renderModeBanner();
    renderSelectedDrugMeta(result.profile.drug);
    renderMechanismExplainCard(result.profile, result.selected_pair);
    renderDiseaseExplainCard(result.profile);
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
  document.getElementById("toast").addEventListener("click", () => {
    const toast = document.getElementById("toast");
    toast.classList.remove("visible");
    if (state.toastTimer) {
      window.clearTimeout(state.toastTimer);
      state.toastTimer = null;
    }
  });

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

  document.getElementById("diseaseSelect").addEventListener("change", async (event) => {
    state.selectedDisease = event.target.value;
    syncCommunityModeUi();
    renderModeBanner();
    if (isCustomMode() || state.selectedDrug) {
      await runStep3Simulation();
      renderScenarioGrid([]);
      await loadPredictions();
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
    state.selectedDisease = "";
    populateDrugSelect();
    renderMicrobePanelNote();
    populateScenarioSelect();
    populateDiseaseSelect();
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
