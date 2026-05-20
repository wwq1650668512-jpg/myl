const state = {
  bootstrap: null,
  mode: "library",
  filteredDrugs: [],
  currentSessionId: null,
};

function $(id) {
  return document.getElementById(id);
}

function formatNumber(value, digits = 2) {
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric.toFixed(digits) : "N/A";
}

function formatPercent(value, digits = 1) {
  const numeric = Number(value);
  return Number.isFinite(numeric) ? `${(numeric * 100).toFixed(digits)}%` : "N/A";
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

async function fetchJson(url, options = {}) {
  const response = await fetch(url, options);
  const payload = await response.json();
  if (!response.ok) {
    throw new Error(payload.error || "Request failed");
  }
  return payload;
}

function setStatus(message) {
  $("statusText").textContent = message;
}

function createStatCard(label, value, note = "") {
  return `
    <div class="stat-card">
      <div class="stat-label">${escapeHtml(label)}</div>
      <div class="stat-value">${escapeHtml(value)}</div>
      <div class="stat-note">${escapeHtml(note)}</div>
    </div>
  `;
}

function createChip(label, value) {
  return `<div class="summary-chip"><strong>${escapeHtml(label)}</strong>${escapeHtml(value)}</div>`;
}

function createTag(label, value) {
  return `<div class="tag"><strong>${escapeHtml(label)}</strong>${escapeHtml(String(value))}</div>`;
}

function createList(items, valueFormatter) {
  if (!items || !items.length) {
    return `<div class="empty-state">暂无结果</div>`;
  }
  return items
    .map((item) => {
      const title = item.microbe_label || item.species_label || item.nt_code || "Unknown";
      const subtitleParts = [];
      if (item.predicted_effect_label) {
        subtitleParts.push(item.predicted_effect_label);
      }
      if (item.predicted_metabolism_label) {
        subtitleParts.push(item.predicted_metabolism_label);
      }
      if (item.nt_code) {
        subtitleParts.push(item.nt_code);
      }
      if (item.delta_abundance !== undefined && item.delta_abundance !== null) {
        subtitleParts.push(`delta ${formatNumber(item.delta_abundance, 4)}`);
      }
      return `
        <div class="list-item">
          <div class="list-item-main">
            <div class="list-item-title">${escapeHtml(title)}</div>
            <div class="list-item-subtitle">${escapeHtml(subtitleParts.join(" · "))}</div>
          </div>
          <div class="list-item-value">${escapeHtml(valueFormatter(item))}</div>
        </div>
      `;
    })
    .join("");
}

function populateScenarios(scenarios) {
  $("scenarioSelect").innerHTML = (scenarios || [])
    .map(
      (item) =>
        `<option value="${escapeHtml(item.scenario_name)}">${escapeHtml(item.scenario_name)}</option>`,
    )
    .join("");
}

function renderDrugOptions() {
  const select = $("drugSelect");
  select.innerHTML = state.filteredDrugs
    .map(
      (drug) =>
        `<option value="${escapeHtml(drug.prestwick_id)}">${escapeHtml(drug.label || drug.chemical_name)}</option>`,
    )
    .join("");
}

function applyDrugFilter() {
  const keyword = $("drugFilterInput").value.trim().toLowerCase();
  const drugs = state.bootstrap?.drugs || [];
  state.filteredDrugs = drugs.filter((drug) => {
    if (!keyword) {
      return true;
    }
    const haystack = `${drug.chemical_name || ""} ${drug.prestwick_id || ""}`.toLowerCase();
    return haystack.includes(keyword);
  });
  renderDrugOptions();
}

function setMode(mode) {
  state.mode = mode;
  $("libraryModeButton").classList.toggle("active", mode === "library");
  $("customModeButton").classList.toggle("active", mode === "custom");
  $("libraryFields").classList.toggle("hidden", mode !== "library");
  $("customFields").classList.toggle("hidden", mode !== "custom");
}

function renderDrugSummary(drug) {
  if (!drug) {
    $("drugSummary").innerHTML = `<div class="empty-state">运行预测后显示药物摘要</div>`;
    return;
  }
  const chips = [
    createChip("名称", drug.chemical_name || "N/A"),
    createChip("ID", drug.prestwick_id || "N/A"),
    createChip("类别", drug.therapeutic_class || "N/A"),
    createChip("作用", drug.therapeutic_effect || "N/A"),
  ];
  if (drug.smiles || drug.canonical_smiles_rdkit) {
    chips.push(createChip("SMILES", drug.canonical_smiles_rdkit || drug.smiles));
  }
  $("drugSummary").innerHTML = chips.join("");
}

function renderStep1(profile) {
  const aggregated = profile?.aggregated || {};
  $("step1Stats").innerHTML = [
    createStatCard("平均 Effect Score", formatNumber(aggregated.mean_predicted_effect_score)),
    createStatCard("平均 Inhibit Prob", formatPercent(aggregated.mean_predicted_inhibit_probability)),
    createStatCard("平均 Promote Prob", formatPercent(aggregated.mean_predicted_promote_probability_refined)),
    createStatCard("置信度", formatNumber(profile?.confidence_score), profile?.confidence_tier || ""),
  ].join("");

  const counts = aggregated.step1_counts || {};
  $("step1Counts").innerHTML = Object.keys(counts).length
    ? Object.entries(counts)
        .map(([label, count]) => createTag(label, count))
        .join("")
    : `<div class="empty-state">暂无统计</div>`;

  $("step1TopList").innerHTML = createList((profile?.top_effect_microbes || []).slice(0, 6), (item) =>
    formatNumber(item.predicted_effect_score),
  );
}

function renderStep2(profile) {
  const aggregated = profile?.aggregated || {};
  $("step2Stats").innerHTML = [
    createStatCard("平均 Metabolized Prob", formatPercent(aggregated.mean_predicted_metabolized_probability)),
    createStatCard("Applicability Rate", formatPercent(aggregated.applicability_rate)),
    createStatCard("平均 Enzyme Support", formatNumber(aggregated.mean_enzyme_support_score)),
    createStatCard("反应类投影对数", String(aggregated.reaction_projection_pairs ?? "N/A")),
  ].join("");

  const counts = aggregated.step2_counts || {};
  $("step2Counts").innerHTML = Object.keys(counts).length
    ? Object.entries(counts)
        .map(([label, count]) => createTag(label, count))
        .join("")
    : `<div class="empty-state">暂无统计</div>`;

  $("step2TopList").innerHTML = createList((profile?.top_metabolism_microbes || []).slice(0, 6), (item) =>
    formatPercent(item.predicted_metabolized_probability),
  );
}

function renderStep3(result) {
  const summary = result?.summary || {};
  $("step3Stats").innerHTML = [
    createStatCard("Development Score", formatNumber(summary.development_score)),
    createStatCard("Health Index", formatNumber(summary.final_health_index)),
    createStatCard("Parent Retention", formatPercent(summary.final_parent_retention_ratio)),
    createStatCard("Metabolite Burden", formatNumber(summary.metabolite_burden_penalty_final)),
  ].join("");

  $("step3Meta").innerHTML = [
    createTag("Scenario", summary.scenario_name || "N/A"),
    createTag("Steps", summary.n_steps ?? "N/A"),
    createTag("Benefit", formatNumber(summary.benefit_subscore_final)),
    createTag("Risk", formatNumber(summary.risk_subscore_final)),
    createTag("Community", summary.community_source || "N/A"),
  ].join("");

  $("step3TopList").innerHTML = createList((result?.top_microbe_changes || []).slice(0, 6), (item) =>
    formatNumber(item.fold_change),
  );
}

function getStep3Payload() {
  return {
    scenario: $("scenarioSelect").value || "healthy_reference",
    n_steps: Number($("nStepsInput").value || 14),
    initial_dose: Number($("initialDoseInput").value || 1.0),
    repeat_dose: Number($("repeatDoseInput").value || 1.0),
  };
}

async function runLibraryPrediction() {
  const drug = $("drugSelect").value;
  if (!drug) {
    throw new Error("请先选择药物。");
  }
  const profile = await fetchJson(`/api/drug-profile?drug=${encodeURIComponent(drug)}`);
  const step3 = await fetchJson("/api/step3/simulate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      drug,
      ...getStep3Payload(),
    }),
  });
  state.currentSessionId = null;
  renderDrugSummary(profile.drug);
  renderStep1(profile);
  renderStep2(profile);
  renderStep3(step3);
}

async function runCustomPrediction() {
  const smiles = $("customSmilesInput").value.trim();
  if (!smiles) {
    throw new Error("请先输入 SMILES。");
  }
  const predictPayload = {
    smiles,
    drug_name: $("customDrugNameInput").value.trim() || undefined,
  };
  const prediction = await fetchJson("/api/custom-drug/predict", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(predictPayload),
  });
  state.currentSessionId = prediction.session_id || null;
  const step3 = await fetchJson("/api/custom-drug/step3/simulate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      session_id: prediction.session_id,
      ...getStep3Payload(),
    }),
  });
  renderDrugSummary(prediction.profile?.drug);
  renderStep1(prediction.profile);
  renderStep2(prediction.profile);
  renderStep3(step3);
}

async function runPrediction() {
  const button = $("runButton");
  button.disabled = true;
  setStatus("正在运行，请稍候…");
  try {
    if (state.mode === "library") {
      await runLibraryPrediction();
    } else {
      await runCustomPrediction();
    }
    setStatus("预测完成。");
  } catch (error) {
    setStatus(`运行失败：${error.message}`);
  } finally {
    button.disabled = false;
  }
}

async function init() {
  try {
    const bootstrap = await fetchJson("/api/bootstrap");
    state.bootstrap = bootstrap;
    state.filteredDrugs = [...(bootstrap.drugs || [])];
    renderDrugOptions();
    populateScenarios(bootstrap.scenarios || []);
    setMode("library");
    setStatus("可以开始运行三步预测。");
    renderDrugSummary(null);
    renderStep1(null);
    renderStep2(null);
    renderStep3(null);
  } catch (error) {
    setStatus(`初始化失败：${error.message}`);
  }
}

$("libraryModeButton").addEventListener("click", () => setMode("library"));
$("customModeButton").addEventListener("click", () => setMode("custom"));
$("drugFilterInput").addEventListener("input", applyDrugFilter);
$("runButton").addEventListener("click", runPrediction);

init();
