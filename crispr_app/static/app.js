let cachedGuides = [];

const toTable = (elementId, rows) => {
  const table = document.getElementById(elementId);
  if (!rows.length) {
    table.innerHTML = "";
    return;
  }
  const cols = Object.keys(rows[0]);
  table.innerHTML = `<thead><tr>${cols.map((c) => `<th>${c}</th>`).join("")}</tr></thead>
    <tbody>${rows.map((r) => `<tr>${cols.map((c) => `<td>${r[c]}</td>`).join("")}</tr>`).join("")}</tbody>`;
};

document.getElementById("designBtn").onclick = async () => {
  const payload = {
    dna_sequence: document.getElementById("dna").value,
    pam: document.getElementById("pam").value,
    guide_length: Number(document.getElementById("guideLength").value),
    min_gc: Number(document.getElementById("minGc").value),
    max_gc: Number(document.getElementById("maxGc").value),
  };

  const resp = await fetch("/api/design", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const data = await resp.json();
  if (!resp.ok) {
    document.getElementById("guideSummary").innerText = data.detail || "Design failed";
    return;
  }

  cachedGuides = data.top_guides;
  document.getElementById("guideSummary").innerText = `Designed ${data.count} guides. Showing top ${Math.min(data.top_guides.length, 100)}.`;
  toTable("guidesTable", data.top_guides.slice(0, 20));
};

document.getElementById("offBtn").onclick = async () => {
  if (!cachedGuides.length) {
    document.getElementById("offSummary").innerText = "Design guides first.";
    return;
  }

  const payload = {
    guides: cachedGuides.map((x) => x.gRNA),
    background_sequence: document.getElementById("background").value,
    max_mismatches: Number(document.getElementById("maxMm").value),
    pam: document.getElementById("pam").value,
  };

  const resp = await fetch("/api/offtargets", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const data = await resp.json();
  if (!resp.ok) {
    document.getElementById("offSummary").innerText = data.detail || "Off-target scan failed";
    return;
  }

  document.getElementById("offSummary").innerText = `Found ${data.count} candidate off-target loci.`;
  toTable("offTable", data.off_targets.slice(0, 30));
};
