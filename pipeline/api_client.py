"""
OMTX Hub API client — using the official omtx Python SDK.

Docs: https://www.omtx.ai/docs/api/hub/routes
SDK:  pip install omtx
Auth: x-api-key header (set OMTX_API_KEY env var)

Pricing (per run, 2026):
  alphafold    $3.00   multimer structure prediction
  bindcraft   $10.00   de novo binder design from receptor PDB    ← Track A leaderboard
  boltz2       $0.85   protein-ligand structure prediction
  boltzgen     $5.00   generative binder/protein design           ← Track A leaderboard
  chai1        $0.85   multimodal structure prediction (CHEAPER than alphafold for screening)
  diffdock     $1.05   diffusion-based ligand docking             ← Track B leaderboard
  flowdock     $0.85   flow-based ligand docking                  ← Track B leaderboard
  neuralplexer $5.00   protein-ligand complex prediction          ← Track B leaderboard
  openfold3    $1.50   multimer + ligand structure prediction
  proteinttt   $1.25   ESMFold test-time tuning
  rfd3         $1.00   all-atom diffusion design
  rosettafold3 $0.75   sequence-to-structure prediction
"""

import os
import time
import uuid
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from omtx import OmClient

# ── Pricing table (USD per run) ───────────────────────────────────────────────
PRICES: dict[str, float] = {
    "alphafold":    3.00,
    "bindcraft":   10.00,
    "boltz2":       0.85,
    "boltzgen":     5.00,
    "chai1":        0.85,
    "diffdock":     1.05,
    "flowdock":     0.85,
    "neuralplexer": 5.00,
    "openfold3":    1.50,
    "proteinttt":   1.25,
    "rfd3":         1.00,
    "rosettafold3": 0.75,
}

# Models that post scores to Hub → Leaderboard
LEADERBOARD_MODELS: dict[str, str] = {
    "bindcraft":   "Track A — de novo binder design",
    "boltzgen":    "Track A — generative binder design",
    "boltz2":      "Track B — protein-ligand ΔG",
    "diffdock":    "Track B — docking confidence",
    "flowdock":    "Track B — docking confidence",
    "neuralplexer":"Track B — co-structure confidence",
}


class JobStatus(str, Enum):
    PENDING   = "pending"
    QUEUED    = "queued"
    RUNNING   = "running"
    COMPLETED = "completed"
    FAILED    = "failed"
    UNKNOWN   = "unknown"

    @classmethod
    def from_str(cls, s: str) -> "JobStatus":
        s = (s or "").lower()
        # OMTX API returns "succeeded" — normalise to "completed"
        if s == "succeeded":
            return cls.COMPLETED
        try:
            return cls(s)
        except ValueError:
            return cls.UNKNOWN


@dataclass
class JobResult:
    job_id: str
    variant_name: str
    model: str
    status: JobStatus
    score: float | None = None
    score_label: str | None = None
    error: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)

    def is_done(self) -> bool:
        return self.status in (JobStatus.COMPLETED, JobStatus.FAILED)

    def __str__(self) -> str:
        score_str = f"{self.score_label}={self.score:.4f}" if self.score is not None else "no score yet"
        return f"[{self.job_id}] {self.variant_name} ({self.model}) — {self.status.value} — {score_str}"


class OmtxClient:
    """
    Wrapper around the official omtx SDK for Track A hackathon submissions.
    Uses v2 Hub routes and the x-api-key header.
    """

    def __init__(self, api_key: str | None = None):
        self.api_key = api_key or os.environ.get("OMTX_API_KEY", "")
        if not self.api_key:
            raise ValueError(
                "OMTX API key required.\n"
                "  export OMTX_API_KEY=your_key_here\n"
                "  Get your key: Hub → User Menu → API Keys"
            )
        self._client = OmClient(api_key=self.api_key)

    # ── Artifact upload ────────────────────────────────────────────────────

    def upload_pdb(self, pdb_path: Path) -> str:
        """Upload a PDB file and return its artifact_id."""
        result = self._client.artifacts.upload(str(pdb_path))
        return result["artifact_id"]

    def upload_cif(self, cif_path: Path) -> str:
        """Upload a CIF file via signed URL (for larger files) and return artifact_id."""
        result = self._client.artifacts.upload_via_signed_url(str(cif_path))
        return result["artifact_id"]

    # ── AlphaFold multimer (receptor + peptide binder complex scoring) ────

    def submit_alphafold(
        self,
        receptor_sequence: str,
        binder_sequence: str,
        variant_name: str,
    ) -> JobResult:
        """
        AlphaFold multimer — $3.00/run.
        receptor_sequence → protein_sequence
        binder_sequence   → protein_sequences_additional
        NOTE: Use submit_chai1() for the same job at $0.85/run.
        """
        idem_key = f"alphafold-{variant_name}-{uuid.uuid4().hex[:8]}"
        data = self._client.hub.submit(
            job_type="hub.alphafold",
            payload={
                "protein_sequence": receptor_sequence,
                "protein_sequences_additional": binder_sequence,
                "model_type": "multimer",
            },
            idempotency_key=idem_key,
        )
        return JobResult(
            job_id=data.get("job_id", "unknown"),
            variant_name=variant_name,
            model="alphafold",
            status=JobStatus.from_str(data.get("status", "queued")),
            raw=data,
        )

    # ── Chai-1 (cheap receptor + peptide screening, $0.85/run) ───────────

    def submit_chai1(
        self,
        receptor_sequence: str,
        binder_sequence: str,
        variant_name: str,
    ) -> JobResult:
        """
        Chai-1 multimodal structure prediction — $0.85/run.
        Same use case as AlphaFold multimer but 3.5x cheaper.
        receptor_sequence → protein_sequence (chain A)
        binder_sequence   → protein_sequences_additional (chain B)
        """
        idem_key = f"chai1-{variant_name}-{uuid.uuid4().hex[:8]}"
        data = self._client.hub.submit(
            job_type="hub.chai1",
            payload={
                "protein_sequence": receptor_sequence,
                "protein_sequences_additional": binder_sequence,
            },
            idempotency_key=idem_key,
        )
        return JobResult(
            job_id=data.get("job_id", "unknown"),
            variant_name=variant_name,
            model="chai1",
            status=JobStatus.from_str(data.get("status", "queued")),
            raw=data,
        )

    # ── Boltz-2 (single-chain folding / ligand) ───────────────────────────

    def submit_boltz2(
        self,
        protein_sequence: str,
        variant_name: str,
        ligand_smiles: str | None = None,
    ) -> JobResult:
        """
        Submit a Boltz-2 job. Boltz-2 Hub API accepts a single protein_sequence
        plus optional ligand_smiles (Track B use case).
        For Track A receptor+peptide complex use submit_alphafold() instead.
        """
        idem_key = f"boltz2-{variant_name}-{uuid.uuid4().hex[:8]}"
        payload: dict[str, Any] = {"protein_sequence": protein_sequence}
        if ligand_smiles:
            payload["ligand_smiles"] = ligand_smiles
        data = self._client.hub.submit(
            job_type="hub.boltz2",
            payload=payload,
            idempotency_key=idem_key,
        )
        return JobResult(
            job_id=data.get("job_id", "unknown"),
            variant_name=variant_name,
            model="boltz2",
            status=JobStatus.from_str(data.get("status", "queued")),
            raw=data,
        )

    # ── BindCraft ──────────────────────────────────────────────────────────

    def submit_bindcraft(
        self,
        receptor_pdb_path: Path,
        variant_name: str,
        target_chains: str = "R",     # chain R matches 7KI0_GLP1R_chainR.pdb
        hotspot_residues: str = "R117,R271,R325,R341",  # chain R prefix
        artifact_id: str | None = None,
    ) -> JobResult:
        """
        Submit a BindCraft job from an uploaded receptor PDB.
        Uploads the PDB if artifact_id is not already cached.
        """
        if artifact_id is None:
            artifact_id = self.upload_pdb(receptor_pdb_path)

        idem_key = f"bindcraft-{variant_name}-{uuid.uuid4().hex[:8]}"
        data = self._client.hub.submit(
            job_type="hub.bindcraft",
            payload={
                "target_pdb_artifact_id": artifact_id,
                "target_chains": target_chains,
                "target_name": f"GLP1R_{variant_name}",
                "hotspot_residues": hotspot_residues,
            },
            idempotency_key=idem_key,
        )
        return JobResult(
            job_id=data.get("job_id", "unknown"),
            variant_name=variant_name,
            model="bindcraft",
            status=JobStatus.from_str(data.get("status", "queued")),
            raw=data,
        )

    # ── BoltzGen ───────────────────────────────────────────────────────────

    def submit_boltzgen(
        self,
        target_cif_artifact_id: str,
        variant_name: str,
        binder_length_min: int = 25,
        binder_length_max: int = 40,
        target_chain_id: str = "F",    # chain F = GLP-1R (entity 6) in the original 7ki0.cif
    ) -> JobResult:
        """Submit a BoltzGen de novo binder design job.

        Requires the *full* RCSB mmCIF (7ki0.cif) — stripped per-chain exports
        lack the _entity_poly_seq table BoltzGen validates.  In 7ki0.cif the
        GLP-1R is label_asym_id F (entity 6).
        """
        idem_key = f"boltzgen-{variant_name}-{uuid.uuid4().hex[:8]}"
        data = self._client.hub.submit(
            job_type="hub.boltzgen",
            payload={
                "protocol": "protein_anything",
                "target_cif_artifact_id": target_cif_artifact_id,
                "target_chain_id": target_chain_id,
                "binder_length_min": binder_length_min,
                "binder_length_max": binder_length_max,
            },
            idempotency_key=idem_key,
        )
        return JobResult(
            job_id=data.get("job_id", "unknown"),
            variant_name=variant_name,
            model="boltzgen",
            status=JobStatus.from_str(data.get("status", "queued")),
            raw=data,
        )

    # ── Job polling ────────────────────────────────────────────────────────

    def get_job(self, job_id: str, variant_name: str = "", model: str = "") -> JobResult:
        """Fetch current status of a job."""
        data = self._client.jobs.status(job_id)
        status = JobStatus.from_str(data.get("status", "unknown"))
        score, label = _extract_score(data)
        return JobResult(
            job_id=job_id,
            variant_name=variant_name or data.get("metadata", {}).get("variant", ""),
            model=model or data.get("job_type", ""),
            status=status,
            score=score,
            score_label=label,
            error=data.get("error"),
            raw=data,
        )

    def poll_until_done(
        self,
        job: JobResult,
        poll_interval: int = 15,
        timeout: int = 1800,
    ) -> JobResult:
        """Block until the job completes or timeout elapses."""
        try:
            data = self._client.jobs.wait(
                job.job_id,
                poll_interval=poll_interval,
                timeout=timeout,
            )
            status = JobStatus.from_str(data.get("status", "unknown"))
            score, label = _extract_score(data)
            job.status = status
            job.score = score
            job.score_label = label
            job.error = data.get("error")
            job.raw = data
        except Exception as e:
            job.status = JobStatus.UNKNOWN
            job.error = str(e)
        return job

    def list_jobs(self, limit: int = 100) -> list[dict[str, Any]]:
        """List recent jobs."""
        return self._client.jobs.history(limit=limit)

    def check_credits(self) -> float:
        """Return available wallet credits (in cents)."""
        profile = self._client.users.profile()
        return profile.get("available_credits", 0)


# ── Score extraction ──────────────────────────────────────────────────────────

def _extract_score(data: dict[str, Any]) -> tuple[float | None, str | None]:
    """
    Try common score paths in OMTX job result payloads.
    Returns (score, label) or (None, None).

    OMTX models return scores in different places:
      - AlphaFold:  response_payload.metrics.iptm_score / ptm_score / plddt_score
                    response_payload.candidates[0].iptm
      - BoltzGen:   response_payload.candidates[0].metrics.iptm
      - BindCraft:  response_payload.candidates[0].metrics.confidence
    Priority: iptm > ptm > confidence > plddt
    """
    resp = data.get("response_payload") or {}

    # 1. Flat metrics dict (AlphaFold / BindCraft style)
    metrics = resp.get("metrics") or {}
    for key, label in [
        ("iptm_score", "iptm"),   # AlphaFold
        ("top_iptm",   "iptm"),   # BindCraft
        ("ptm_score",  "ptm"),    # AlphaFold
        ("top_plddt",  "plddt"),  # BindCraft fallback
        ("plddt_score","plddt"),  # AlphaFold fallback
        ("confidence", "confidence"),
    ]:
        val = metrics.get(key)
        if val is not None:
            try:
                return float(val), label
            except (TypeError, ValueError):
                pass

    # 2. Candidates list (BoltzGen / BindCraft / AlphaFold candidates array)
    candidates = resp.get("candidates") or []
    if candidates:
        c = candidates[0]
        cmetrics = c if isinstance(c, dict) else {}
        # BoltzGen nests under metrics sub-dict; AlphaFold is flat
        inner = cmetrics.get("metrics") or cmetrics
        for key in ("iptm", "ptm", "confidence", "score", "plddt"):
            val = inner.get(key)
            if val is not None:
                try:
                    return float(val), key
                except (TypeError, ValueError):
                    pass

    # 3. Fallback: scan top-level and result/results sub-dicts
    results = data.get("result") or data.get("results") or {}
    for d in [data, results, resp, metrics]:
        for key in ("confidence", "score", "plddt", "ipae", "ptm",
                    "iptm", "dg", "delta_g", "binding_energy"):
            val = d.get(key)
            if val is not None:
                try:
                    return float(val), key
                except (TypeError, ValueError):
                    pass
    return None, None
