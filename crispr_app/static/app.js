"use strict";

const $ = (id) => document.getElementById(id);
const EXAMPLE_SEQ =
  "ATGGCCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCC" +
  "GCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAG" +
  "CTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCCGCG" +
  "GTGGCGGTCTGGACCACGCCGGAGAGCGTCGAAGCGGGGGCGGTGTTCGCCGAGATCGGCCCGCGCATGGCC";

let lastGuides = [];
let lastOff = [];
let lastSpec = [];
let lastPeg = [];

/* ---------- helpers ---------- */

function badge(value) {
  const v = Number(value);
  const cls = v >= 0.6 ? "good" : v >= 0.4 ? "mid" : "low";
  return `<span class="badge ${cls}">${v.toFixed(3)}</span>`;
}

function specBadgeClass(v) {
  return v >= 80 ? "good" : v >= 50 ? "mid" : "low";
}

function seqWithPam(guide, pam) {
  // Render guide in mono with the PAM appended and highlighted.
  return `<span class="seq">${guide}<span class="pam">${pam || ""}</span></span>`;
}

function copyBtn(text) {
  return `<button class="iconbtn" title="Copy" data-copy="${text}">⧉</button>`;
}

async function postJSON(url, body) {
  const resp = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const data = await resp.json().catch(() => ({}));
  if (!resp.ok) throw new Error(data.detail || `Request failed (${resp.status})`);
  return data;
}

function busy(btn, on, label) {
  if (on) {
    btn.dataset.label = btn.innerHTML;
    btn.disabled = true;
    btn.innerHTML = `<span class="spinner"></span>${label || "Working…"}`;
  } else {
    btn.disabled = false;
    btn.innerHTML = btn.dataset.label || btn.innerHTML;
  }
}

function setStatus(el, msg, isError = false) {
  el.textContent = msg;
  el.classList.toggle("error", isError);
}

function toCSV(rows) {
  if (!rows.length) return "";
  const cols = Object.keys(rows[0]);
  const esc = (v) => `"${String(v).replace(/"/g, '""')}"`;
  return [cols.join(","), ...rows.map((r) => cols.map((c) => esc(r[c])).join(","))].join("\n");
}

function toTable(elementId, rows) {
  const t = $(elementId);
  if (!rows || !rows.length) { t.innerHTML = ""; return; }
  const cols = Object.keys(rows[0]);
  t.innerHTML =
    `<thead><tr>${cols.map((c) => `<th>${c}</th>`).join("")}</tr></thead>` +
    `<tbody>${rows.map((r) => `<tr>${cols.map((c) => `<td>${r[c]}</td>`).join("")}</tr>`).join("")}</tbody>`;
}

function download(filename, text) {
  const blob = new Blob([text], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}

/* ---------- rendering ---------- */

function renderGuides(rows) {
  const t = $("guidesTable");
  if (!rows.length) { t.innerHTML = `<tbody><tr><td class="empty">No guides matched the GC / quality filters.</td></tr></tbody>`; return; }
  const head = `<thead><tr>
    <th>#</th><th>Strand</th><th>Start</th><th>Guide (5'→3') + PAM</th>
    <th>GC%</th><th title="0–100 prioritization score combining on-target predictors. Relative ranking, not a literal % editing rate.">Score</th><th></th></tr></thead>`;
  const score100 = (v) => {
    const n = Math.round(Number(v) * 100);
    const cls = n >= 60 ? "good" : n >= 40 ? "mid" : "low";
    return `<span class="badge ${cls}">${n}</span>`;
  };
  const body = rows.map((r, i) => `<tr>
    <td>${i + 1}</td>
    <td>${r.Strand}</td>
    <td>${r.Start}</td>
    <td>${seqWithPam(r.gRNA, r.PAM)} ${copyBtn(r.gRNA)}</td>
    <td>${r["GC%"]}</td>
    <td>${score100(r.ConsensusScore)}</td>
    <td><button class="iconbtn" data-explain="${r.gRNA}" title="Show the score breakdown">Details</button></td>
  </tr>`).join("");
  t.innerHTML = head + `<tbody>${body}</tbody>`;
}

function renderSpecificity(spec) {
  const g = $("specGrid");
  if (!spec || !spec.length) { g.innerHTML = ""; return; }
  g.innerHTML = spec.map((s) => `
    <div class="spec">
      <div class="g">${s.gRNA}</div>
      <div class="nums">
        <div>CFD spec<b class="badge ${specBadgeClass(s.CFD_Specificity)}">${s.CFD_Specificity}</b></div>
        <div>MIT spec<b class="badge ${specBadgeClass(s.MIT_Specificity)}">${s.MIT_Specificity}</b></div>
        <div>Off-targets<b>${s.OffTargetCount}</b></div>
      </div>
    </div>`).join("");
}

function renderOff(rows) {
  const t = $("offTable");
  if (!rows.length) { t.innerHTML = `<tbody><tr><td class="empty">No off-target loci within the mismatch threshold — a good sign.</td></tr></tbody>`; return; }
  const head = `<thead><tr><th>Guide</th><th>Strand</th><th>Pos</th><th>MM</th><th>Off-target + PAM</th><th>CFD</th><th>MIT</th></tr></thead>`;
  const body = rows.slice(0, 50).map((r) => `<tr>
    <td class="seq">${r.gRNA}</td>
    <td>${r.Strand}</td>
    <td>${r.OffTargetPos}</td>
    <td>${r.Mismatches}</td>
    <td>${seqWithPam(r.TargetSeq, r.PAM)}</td>
    <td>${badge(r.CFD_Score)}</td>
    <td>${Number(r.MIT_Score).toFixed(1)}</td>
  </tr>`).join("");
  t.innerHTML = head + `<tbody>${body}</tbody>`;
}

function renderPeg(rows) {
  const t = $("primeTable");
  if (!rows.length) { t.innerHTML = `<tbody><tr><td class="empty">No pegRNA found — try a target index nearer an NGG PAM.</td></tr></tbody>`; return; }
  const head = `<thead><tr><th>#</th><th>Spacer</th><th>PAM</th><th>Nick</th><th>PBS</th><th>PBS Tm</th><th>RTT</th><th>Score</th></tr></thead>`;
  const body = rows.slice(0, 30).map((r, i) => `<tr>
    <td>${i + 1}</td>
    <td class="seq">${r.Spacer} ${copyBtn(r.Spacer)}</td>
    <td>${r.PAM}</td>
    <td>${r.NickPos}</td>
    <td class="seq">${r.PBS_Seq}</td>
    <td>${r.PBS_Tm}</td>
    <td class="seq">${r.RTT_Seq}</td>
    <td>${badge(r.Score)}</td>
  </tr>`).join("");
  t.innerHTML = head + `<tbody>${body}</tbody>`;
}

function renderExplain(data) {
  const p = $("explainPanel");
  const c = data.contributions || {};
  const entries = Object.entries(c);
  const max = Math.max(0.3, ...entries.map(([, v]) => Math.abs(v)));
  const bars = entries.map(([name, v]) => {
    const pct = (Math.abs(v) / max) * 50;
    const side = v >= 0 ? `left:50%;width:${pct}%` : `right:50%;width:${pct}%`;
    return `<div class="bar-row">
      <span class="name">${name.replace(/_/g, " ")}</span>
      <span class="bar-track"><span class="bar-mid"></span><span class="bar-fill ${v >= 0 ? "pos" : "neg"}" style="${side}"></span></span>
      <span class="val">${v >= 0 ? "+" : ""}${v.toFixed(3)}</span>
    </div>`;
  }).join("");
  let ciHtml = "";
  if (data.interval) {
    const ci = data.interval;
    ciHtml = `<p class="status"><b>${Math.round(ci.coverage * 100)}% confidence interval</b> (conformal, guaranteed coverage):
      on-target ${ci.point} &nbsp;[${ci.low} – ${ci.high}]. Wide intervals = genuine model uncertainty.</p>`;
  }
  p.innerHTML = `<h4>Why ${data.guide} scored ${data.score} &nbsp;
      <span class="status" style="display:inline">(GC ${data.gc_percent}% · Tm ${data.tm_celsius} °C)</span></h4>
    <div class="bars">${bars}</div>
    ${ciHtml}
    <p class="status">Bars are additive log-odds contributions; the sum is squashed to a 0–1 efficiency probability.</p>`;
  p.classList.add("open");
  p.scrollIntoView({ behavior: "smooth", block: "nearest" });
}

/* ---------- actions ---------- */

(async () => {
  const friendly = { linear: "trained model", heuristic: "built-in (interpretable)", onnx: "deep model", crisprscan: "CRISPRscan" };
  try {
    const m = await (await fetch("/api/models")).json();
    $("modelTag").textContent = `on-target model: ${friendly[m.active] || m.active}`;
    $("modelTag").title = `available: ${m.available.map((x) => friendly[x] || x).join(", ")}`;
  } catch (_) { $("modelTag").textContent = "on-target model: built-in"; }
})();

$("exampleBtn").onclick = () => { $("dna").value = EXAMPLE_SEQ; $("background").value = EXAMPLE_SEQ; };
$("clearBtn").onclick = () => { $("dna").value = ""; };

$("fastaBtn").onclick = async () => {
  const contents = $("dna").value.trim();
  if (!contents) return;
  try {
    const data = await postJSON("/api/upload-fasta", { contents });
    $("dna").value = data.sequence;
    setStatus($("guideSummary"), `Parsed ${data.length} bp${data.message ? " — " + data.message : ""}.`);
  } catch (e) { setStatus($("guideSummary"), e.message, true); }
};

$("designBtn").onclick = async () => {
  const btn = $("designBtn");
  busy(btn, true, "Designing…");
  try {
    const data = await postJSON("/api/design", {
      dna_sequence: $("dna").value,
      pam: $("pam").value,
      guide_length: Number($("guideLength").value),
      min_gc: Number($("minGc").value),
      max_gc: Number($("maxGc").value),
    });
    lastGuides = data.top_guides;
    renderGuides(lastGuides);
    setStatus($("guideSummary"), `Designed ${data.count} candidate guides — showing top ${Math.min(lastGuides.length, 100)}, ranked by consensus score.`);
    $("exportGuidesBtn").disabled = !lastGuides.length;
    $("offBtn").disabled = !lastGuides.length;
    $("explainPanel").classList.remove("open");

    // Populate the edit-outcome guide picker with forward-strand guides
    // (those are present in the pasted sequence, so the cut site is locatable).
    const fwd = lastGuides.filter((g) => g.Strand === "+");
    const sel = $("simGuide");
    sel.innerHTML = fwd.length
      ? fwd.map((g) => `<option value="${g.gRNA}">${g.gRNA} (start ${g.Start}, score ${Math.round(g.ConsensusScore * 100)})</option>`).join("")
      : `<option value="">— no + strand guides —</option>`;
    $("simBtn").disabled = !fwd.length;
  } catch (e) { setStatus($("guideSummary"), e.message, true); }
  finally { busy(btn, false); }
};

$("offBtn").onclick = async () => {
  if (!lastGuides.length) { setStatus($("offSummary"), "Design guides first.", true); return; }
  const btn = $("offBtn");
  busy(btn, true, "Scanning…");
  try {
    const data = await postJSON("/api/offtargets", {
      guides: lastGuides.map((g) => g.gRNA),
      background_sequence: $("background").value,
      max_mismatches: Number($("maxMm").value),
      pam: $("pam").value,
    });
    lastOff = data.off_targets || [];
    lastSpec = data.specificity || [];
    renderSpecificity(lastSpec);
    renderOff(lastOff);
    setStatus($("offSummary"), `${data.count} candidate off-target loci across both strands. Specificity (0–100, higher is better) summarised per guide below.`);
    $("exportOffBtn").disabled = !lastOff.length;
  } catch (e) { setStatus($("offSummary"), e.message, true); }
  finally { busy(btn, false); }
};

$("primeBtn").onclick = async () => {
  const btn = $("primeBtn");
  busy(btn, true, "Designing…");
  try {
    const data = await postJSON("/api/prime-design", {
      dna_sequence: $("dna").value,
      target_pos: Number($("primePos").value),
      desired_edit: $("primeBase").value,
    });
    lastPeg = data.pegRNAs || [];
    renderPeg(lastPeg);
    setStatus($("primeSummary"), `${data.count} pegRNA candidates — showing top ${Math.min(lastPeg.length, 30)} by PRIDICT2.0-informed score.`);
    $("exportPrimeBtn").disabled = !lastPeg.length;
  } catch (e) { setStatus($("primeSummary"), e.message, true); }
  finally { busy(btn, false); }
};

$("simBtn").onclick = async () => {
  const guide = $("simGuide").value;
  if (!guide) { setStatus($("simSummary"), "Design guides first, then pick one.", true); return; }
  const btn = $("simBtn");
  busy(btn, true, "Simulating…");
  try {
    const data = await postJSON("/api/simulate", {
      dna_sequence: $("dna").value,
      guide,
      edit_offset: 17,           // SpCas9 cut ≈ 3 bp 5' of the PAM
      edit_type: $("editType").value,
    });
    const flag = (on, label) =>
      `<span class="badge ${on ? "low" : "good"}">${label}: ${on ? "yes" : "no"}</span>`;
    $("simPanel").innerHTML = `
      <h4>Predicted protein consequence &nbsp; ${flag(data.frameshift, "frameshift")} ${flag(data.stop_lost, "stop lost")}</h4>
      <p class="status">Protein before (to first stop):</p><div class="seq">${data.protein_before || "—"}</div>
      <p class="status" style="margin-top:.6rem">Protein after the edit:</p><div class="seq">${data.protein_after || "—"}</div>`;
    $("simPanel").classList.add("open");
    setStatus($("simSummary"), `Simulated ${$("editType").options[$("editType").selectedIndex].text.toLowerCase()} at the cut site.`);
    toTable("indelTable", data.indel_panel);
  } catch (e) { setStatus($("simSummary"), e.message, true); $("simPanel").classList.remove("open"); }
  finally { busy(btn, false); }
};

/* export */
$("exportGuidesBtn").onclick = () => download("guides.csv", toCSV(lastGuides));
$("exportOffBtn").onclick = () => download("offtargets.csv", toCSV(lastOff));
$("exportPrimeBtn").onclick = () => download("pegrnas.csv", toCSV(lastPeg));

/* delegated: copy + explain */
document.addEventListener("click", async (e) => {
  const copy = e.target.closest("[data-copy]");
  if (copy) {
    try { await navigator.clipboard.writeText(copy.dataset.copy); copy.textContent = "✓"; setTimeout(() => (copy.textContent = "⧉"), 900); } catch (_) {}
    return;
  }
  const ex = e.target.closest("[data-explain]");
  if (ex) {
    try { renderExplain(await postJSON("/api/explain", { guide: ex.dataset.explain })); }
    catch (err) { setStatus($("guideSummary"), err.message, true); }
  }
});
