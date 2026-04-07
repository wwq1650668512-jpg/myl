from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from enum import Enum
from typing import Optional


class EffectLabel(str, Enum):
    INHIBIT = "inhibit"
    PROMOTE = "promote"
    NO_EFFECT = "no_effect"


class MetabolismLabel(str, Enum):
    METABOLIZED = "metabolized"
    NOT_METABOLIZED = "not_metabolized"
    UNCERTAIN = "uncertain"


class ReactionClass(str, Enum):
    REDUCTION = "reduction"
    HYDROLYSIS = "hydrolysis"
    DEACETYLATION = "deacetylation"
    DEHYDROXYLATION = "dehydroxylation"
    DEMETHYLATION = "demethylation"
    DECONJUGATION = "deconjugation"
    RING_CLEAVAGE = "ring_cleavage"
    BIOACCUMULATION_OR_UNRESOLVED_DEPLETION = "bioaccumulation_or_unresolved_depletion"
    OTHER = "other"


@dataclass(slots=True)
class DrugEntity:
    primary_id: str
    name: str
    canonical_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    pubchem_cid: Optional[str] = None
    chebi_id: Optional[str] = None
    aliases: tuple[str, ...] = ()
    properties: dict[str, float | int | str] = field(default_factory=dict)


@dataclass(slots=True)
class MicrobeEntity:
    primary_id: str
    name: str
    ncbi_taxonomy_id: Optional[str] = None
    rank: Optional[str] = None
    strain: Optional[str] = None
    lineage: tuple[str, ...] = ()
    genome_features: dict[str, float | int | str] = field(default_factory=dict)


@dataclass(slots=True)
class DrugMicrobeEffectRecord:
    drug_id: str
    microbe_id: str
    label: EffectLabel
    effect_score: Optional[float] = None
    assay_name: Optional[str] = None
    source_id: Optional[str] = None
    concentration: Optional[float] = None
    concentration_unit: Optional[str] = None
    p_value: Optional[float] = None
    q_value: Optional[float] = None
    metadata: dict[str, float | int | str | bool] = field(default_factory=dict)

    def is_quantitative(self) -> bool:
        return self.effect_score is not None


@dataclass(slots=True)
class MicrobeDrugMetabolismRecord:
    drug_id: str
    microbe_id: str
    label: MetabolismLabel
    reaction_class: Optional[ReactionClass] = None
    source_id: Optional[str] = None
    parent_depletion_fraction: Optional[float] = None
    product_ids: tuple[str, ...] = ()
    evidence_gene_ids: tuple[str, ...] = ()
    metadata: dict[str, float | int | str | bool] = field(default_factory=dict)

    def has_resolved_products(self) -> bool:
        return len(self.product_ids) > 0


@dataclass(slots=True)
class CommunityState:
    sample_id: str
    abundances: dict[str, float]
    parent_drugs: dict[str, float] = field(default_factory=dict)
    metabolites: dict[str, float] = field(default_factory=dict)
    health_index: Optional[float] = None
    metadata: dict[str, float | int | str | bool] = field(default_factory=dict)

    def normalized_abundances(self) -> dict[str, float]:
        total = sum(value for value in self.abundances.values() if value > 0)
        if total <= 0:
            return {}
        return {
            key: value / total
            for key, value in self.abundances.items()
            if value > 0
        }


@dataclass(slots=True)
class SimulationResult:
    initial_state: CommunityState
    final_state: CommunityState
    timepoints: tuple[int, ...] = ()
    state_history: tuple[CommunityState, ...] = ()
    development_score: Optional[float] = None
    uncertainty: dict[str, float] = field(default_factory=dict)
    notes: tuple[str, ...] = ()

