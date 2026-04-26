// ── Drag & drop + file upload ─────────────────────────────────────────────────

const dropZone   = document.getElementById('dropZone');
const fileInput  = document.getElementById('fileInput');
const seqInput   = document.getElementById('sequenceInput');
const fileNameEl = document.getElementById('fileName');

dropZone.addEventListener('dragover',  e => { e.preventDefault(); dropZone.classList.add('over'); });
dropZone.addEventListener('dragleave', ()  => dropZone.classList.remove('over'));
dropZone.addEventListener('drop', e => {
  e.preventDefault();
  dropZone.classList.remove('over');
  handleFile(e.dataTransfer.files[0]);
});
dropZone.addEventListener('click', e => {
  if (e.target === dropZone || e.target.tagName === 'P' || e.target.tagName === 'SPAN') {
    fileInput.click();
  }
});
fileInput.addEventListener('change', () => handleFile(fileInput.files[0]));

function handleFile(file) {
  if (!file) return;
  const reader = new FileReader();
  reader.onload = e => {
    seqInput.value = e.target.result;
    fileNameEl.textContent = `📄 ${file.name} loaded`;
    autoDetectStrand(e.target.result);
  };
  reader.readAsText(file);
}

// Live strand type reveal
seqInput.addEventListener('input', () => autoDetectStrand(seqInput.value));

function autoDetectStrand(raw) {
  const seq = raw.replace(/>/g, '\n').split('\n')
    .filter(l => !l.startsWith('>'))
    .join('').replace(/\s/g,'').toUpperCase();
  const hasT = /T/.test(seq);
  const hasU = /U/.test(seq);
  const strandGroup = document.getElementById('strandGroup');
  // Show strand type choice only for pure DNA (has T, no U)
  if (hasT && !hasU && seq.length > 0) {
    strandGroup.style.display = 'block';
  } else {
    strandGroup.style.display = 'none';
  }
}


// ── Main analyze function ─────────────────────────────────────────────────────

async function analyze() {
  const seq = seqInput.value.trim();
  const errBox = document.getElementById('inputError');
  errBox.style.display = 'none';

  if (!seq) {
    errBox.textContent = '⚠ Please enter or upload a sequence first.';
    errBox.style.display = 'block';
    return;
  }

  const strandType = document.querySelector('input[name="strand"]:checked')?.value || 'non-template';

  // Show loading
  const overlay = document.getElementById('loadingOverlay');
  const loadMsg = document.getElementById('loadingMsg');
  overlay.style.display = 'flex';
  loadMsg.textContent = 'Analyzing sequence…';

  const loadingMsgs = [
    'Detecting sequence type…',
    'Running transcription…',
    'Translating codons…',
    'Characterizing protein…',
    'Querying UniProt database…',
  ];
  let msgIdx = 0;
  const msgInterval = setInterval(() => {
    loadMsg.textContent = loadingMsgs[Math.min(++msgIdx, loadingMsgs.length - 1)];
  }, 1800);

  try {
    const formData = new FormData();
    formData.append('sequence', seq);
    formData.append('strand_type', strandType);

    const resp = await fetch('/analyze', { method: 'POST', body: formData });
    const data = await resp.json();

    clearInterval(msgInterval);
    overlay.style.display = 'none';

    if (data.error) {
      errBox.textContent = `⚠ ${data.error}`;
      errBox.style.display = 'block';
      return;
    }

    renderResults(data);

  } catch (err) {
    clearInterval(msgInterval);
    overlay.style.display = 'none';
    errBox.textContent = `⚠ Network error: ${err.message}`;
    errBox.style.display = 'block';
  }
}


// ── Render all results ────────────────────────────────────────────────────────

function renderResults(data) {
  document.getElementById('results').style.display = 'flex';
  document.getElementById('results').style.flexDirection = 'column';
  document.getElementById('results').style.gap = '24px';

  renderDetection(data);

  if (data.detection.type === 'INVALID') return;

  renderTranscription(data);
  renderTranslation(data);
  renderAmino(data);
  renderProtein(data);

  // Scroll to results
  document.getElementById('results').scrollIntoView({ behavior: 'smooth' });
}


// ── Step 1: Detection ─────────────────────────────────────────────────────────

function renderDetection(data) {
  const d = data.detection;
  const el = document.getElementById('detectionResult');
  document.getElementById('cardDetection').style.display = 'block';

  let badgeClass, badgeText;
  if      (d.type === 'DNA')     { badgeClass = 'badge-dna';     badgeText = '🔵 DNA Detected'; }
  else if (d.type === 'RNA')     { badgeClass = 'badge-rna';     badgeText = '🟡 RNA Detected'; }
  else                            { badgeClass = 'badge-invalid'; badgeText = '🔴 INVALID Sequence'; }

  let html = `
    <div style="margin-bottom:12px;">
      <span class="badge ${badgeClass}" style="font-size:.95rem; padding:6px 16px;">${badgeText}</span>
    </div>
    <span class="seq-label">Input sequence (${data.sequence.length} bases):</span>
    <div class="seq-block">${highlightSeq(data.sequence, d.type)}</div>
    <div class="explanation"><strong>How was this detected?</strong><br/>${d.explanation}</div>
  `;

  if (d.type === 'INVALID' && d.invalid_chars?.length) {
    html += `<div class="error-box">Invalid characters found: <strong>${d.invalid_chars.join(', ')}</strong></div>`;
  }

  el.innerHTML = html;
}


// ── Step 2: Transcription ─────────────────────────────────────────────────────

function renderTranscription(data) {
  const t = data.transcription;
  const el = document.getElementById('transcriptionResult');
  document.getElementById('cardTranscription').style.display = 'block';

  const strandLabel = data.seq_type === 'RNA' ? '' :
    `<p style="margin-bottom:8px; font-size:.85rem; color:#555;">
       Strand type provided: <strong>${data.strand_type === 'non-template' ? 'Non-template (coding) strand' : 'Template (antisense) strand'}</strong>
     </p>`;

  el.innerHTML = `
    ${strandLabel}
    <span class="seq-label">Input sequence:</span>
    <div class="seq-block">${highlightSeq(t.input_seq, data.seq_type)}</div>
    <div style="text-align:center; font-size:1.3rem; color:#1a6b4a; margin:4px 0;">↓</div>
    <span class="seq-label">mRNA sequence (${t.mrna.length} bases):</span>
    <div class="seq-block" style="color:#ffd060;">${highlightRNA(t.mrna)}</div>
    <div class="explanation"><strong>What is transcription?</strong><br/>${t.explanation}</div>
  `;
}


// ── Step 3: Translation ───────────────────────────────────────────────────────

function renderTranslation(data) {
  const tr = data.translation;
  const el  = document.getElementById('translationResult');
  document.getElementById('cardTranslation').style.display = 'block';

  if (!tr.found_start) {
    el.innerHTML = `<div class="explanation">${tr.explanation}</div>`;
    return;
  }

  let rows = '';
  tr.codons.forEach((c, i) => {
    const cls = c.is_start ? 'is-start' : c.is_stop ? 'is-stop' : '';
    const badge = c.is_start ? '<span class="badge badge-start">START</span>' :
                  c.is_stop  ? '<span class="badge badge-stop">STOP</span>'  : '';
    rows += `
      <tr class="${cls}">
        <td>${i + 1}</td>
        <td class="mono">${c.codon}</td>
        <td>${c.name}</td>
        <td class="mono">${c.three}</td>
        <td class="mono">${c.one}</td>
        <td>${badge}</td>
      </tr>`;
  });

  el.innerHTML = `
    <div class="explanation"><strong>What is translation?</strong><br/>${tr.explanation}</div>
    <div class="codon-table-wrapper">
      <table class="codon-table">
        <thead>
          <tr><th>#</th><th>Codon</th><th>Amino Acid</th><th>3-Letter</th><th>1-Letter</th><th>Note</th></tr>
        </thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  `;
}


// ── Step 4: Amino Acids ───────────────────────────────────────────────────────

function renderAmino(data) {
  const aa  = data.translation.amino_acids;
  const el  = document.getElementById('aminoResult');
  document.getElementById('cardAmino').style.display = 'block';

  const chips = aa.map(a => {
    const isStop = a.name === 'STOP';
    return `
      <div class="aa-chip ${isStop ? 'stop-chip' : ''}">
        <div class="aa-one">${a.one}</div>
        <div class="aa-three">${a.three}</div>
        <div class="aa-name">${a.name}</div>
      </div>`;
  }).join('');

  const explanation = `
    <strong>What are amino acids?</strong><br/>
    Amino acids are the building blocks of proteins. Each codon in the mRNA corresponds to a specific amino acid.
    When amino acids are joined together in a chain, the result is called a <strong>polypeptide chain</strong>.
    This chain then folds into a specific 3D shape to become a functional protein.
    The sequence and identity of amino acids determines exactly what the protein does in the body.
    The chain shown here contains <strong>${aa.filter(a => a.name !== 'STOP').length}</strong> amino acid(s).
  `;

  el.innerHTML = `
    <div class="explanation">${explanation}</div>
    <div class="aa-chain">${chips}</div>
    <p style="font-size:.82rem; color:#888; margin-top:8px;">
      Each tile shows the one-letter code (top), three-letter abbreviation (middle), and full name (bottom).
      <span class="badge badge-start">Green</span> = methionine (start) · 
      <span class="badge badge-stop">Red</span> = stop signal
    </p>
  `;
}


// ── Step 5: Protein ───────────────────────────────────────────────────────────

function renderProtein(data) {
  const pc  = data.protein_char;
  const db  = data.db_result;
  const seq = data.translation.protein_seq;
  const el  = document.getElementById('proteinResult');
  document.getElementById('cardProtein').style.display = 'block';

  // Stats
  const stats = [
    { value: pc.length,           label: 'Amino Acids' },
    { value: pc.hydrophobic,      label: 'Hydrophobic' },
    { value: pc.polar,            label: 'Polar' },
    { value: pc.positive_charged, label: 'Pos. Charged (+)' },
    { value: pc.negative_charged, label: 'Neg. Charged (−)' },
  ].map(s => `
    <div class="stat-box">
      <div class="stat-value">${s.value}</div>
      <div class="stat-label">${s.label}</div>
    </div>`).join('');

  // Amino acid composition table
  const compRows = Object.entries(pc.composition).sort((a,b) => b[1]-a[1]).map(([name, count]) =>
    `<tr><td>${name}</td><td>${count}</td></tr>`
  ).join('');

  // DB results
  let dbHtml = '';
  if (db.success && db.results.length) {
    dbHtml = `
      <h3 style="margin: 18px 0 10px; font-size:1rem; color:#0d3b2e;">
        🌐 Database Results — ${db.source}
      </h3>
      <div class="explanation">
        <strong>What do these results mean?</strong><br/>
        The protein sequence was submitted to the UniProt database — one of the world's largest
        collections of protein information. The results below show real proteins that are similar
        to the one produced from your sequence, what organism they come from, and what they do.
        A higher identity score means the match is more closely related to your protein.
      </div>`;
    db.results.forEach((hit, i) => {
      dbHtml += `
        <div class="db-hit">
          <h4>#${i+1} · <a href="${hit.url}" target="_blank">${hit.accession}</a> — ${hit.protein || hit.name}</h4>
          <p>🦠 <strong>Organism:</strong> ${hit.organism}</p>
          <p>🔬 <strong>Identity:</strong> ${hit.identity}</p>
          <p>📋 <strong>Function:</strong> ${hit.function}</p>
          <p class="db-note"><a href="${hit.url}" target="_blank">View full entry on UniProt →</a></p>
        </div>`;
    });
    if (db.note) dbHtml += `<p class="db-note">Note: ${db.note}</p>`;
  } else {
    dbHtml = `
      <h3 style="margin:18px 0 10px; font-size:1rem; color:#0d3b2e;">🌐 Database Lookup</h3>
      <div class="explanation">
        <strong>What do these results mean?</strong><br/>
        The protein sequence was submitted to UniProt for matching against known proteins.
        ${db.message || ''}
        ${db.manual_url ? `<br/>You can manually BLAST your sequence at: <a href="${db.manual_url}" target="_blank">${db.manual_url}</a>` : ''}
      </div>`;
  }

  el.innerHTML = `
    <div class="explanation"><strong>What is a protein?</strong><br/>${pc.explanation}</div>
    ${seq ? `
      <p style="font-size:.85rem; font-weight:700; color:#0d3b2e; margin-bottom:4px;">Protein sequence (one-letter codes):</p>
      <div class="protein-seq">${seq}</div>
    ` : ''}
    <h3 style="margin:16px 0 10px; font-size:1rem; color:#0d3b2e;">📊 Composition &amp; Properties</h3>
    <div class="stats-grid">${stats}</div>
    <p style="font-size:.85rem; color:#555; margin:4px 0;"><strong>Net charge:</strong> ${pc.charge_desc}</p>
    ${compRows ? `
      <table class="composition-table" style="margin-top:12px;">
        <thead><tr><th>Amino Acid</th><th>Count</th></tr></thead>
        <tbody>${compRows}</tbody>
      </table>` : ''}
    ${dbHtml}
  `;
}


// ── Sequence highlighting helpers ─────────────────────────────────────────────

const DNA_COLORS = { A: '#60a5fa', T: '#f87171', C: '#4ade80', G: '#facc15' };
const RNA_COLORS = { A: '#60a5fa', U: '#fb923c', C: '#4ade80', G: '#facc15' };

function highlightSeq(seq, type) {
  const colors = type === 'RNA' ? RNA_COLORS : DNA_COLORS;
  return seq.split('').map(b =>
    colors[b] ? `<span style="color:${colors[b]}">${b}</span>` : b
  ).join('');
}

function highlightRNA(seq) {
  return seq.split('').map(b =>
    RNA_COLORS[b] ? `<span style="color:${RNA_COLORS[b]}">${b}</span>` : b
  ).join('');
}
