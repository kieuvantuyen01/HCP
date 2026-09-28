use std::fmt;
use std::str::FromStr;

/// Named, reproducible feature configurations used by the ablation study.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AblationConfig {
    Full,
    NoDecomposition,
    NoDispatch,
    NoSeries,
    NoTwoCut,
    NoRepair,
    UndirectedCubic,
    OneAlternatingWorker,
}

impl Default for AblationConfig {
    fn default() -> Self {
        Self::Full
    }
}

impl AblationConfig {
    pub const ALL: [Self; 8] = [
        Self::Full,
        Self::NoDecomposition,
        Self::NoDispatch,
        Self::NoSeries,
        Self::NoTwoCut,
        Self::NoRepair,
        Self::UndirectedCubic,
        Self::OneAlternatingWorker,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::NoDecomposition => "no-decomposition",
            Self::NoDispatch => "no-dispatch",
            Self::NoSeries => "no-series",
            Self::NoTwoCut => "no-two-cut",
            Self::NoRepair => "no-repair",
            Self::UndirectedCubic => "undirected-cubic",
            Self::OneAlternatingWorker => "one-alternating-worker",
        }
    }

    pub fn dispatch_enabled(self) -> bool {
        !matches!(self, Self::NoDecomposition | Self::NoDispatch)
    }

    pub fn series_enabled(self) -> bool {
        !matches!(self, Self::NoDecomposition | Self::NoSeries)
    }

    pub fn two_cut_enabled(self) -> bool {
        !matches!(self, Self::NoDecomposition | Self::NoTwoCut)
    }

    pub fn repair_enabled(self) -> bool {
        !matches!(self, Self::NoRepair)
    }

    pub fn directed_cubic_enabled(self) -> bool {
        self.dispatch_enabled() && !matches!(self, Self::UndirectedCubic)
    }

    /// The alternating frontend contracts forced series chains internally, so it
    /// cannot be selected in the `no-series` condition.
    pub fn alternating_dispatch_enabled(self) -> bool {
        self.dispatch_enabled() && self.series_enabled()
    }

    pub fn alternating_workers(self) -> usize {
        if matches!(self, Self::OneAlternatingWorker) {
            1
        } else {
            3
        }
    }
}

impl fmt::Display for AblationConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for AblationConfig {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let normalized = value.trim().to_ascii_lowercase().replace('_', "-");
        match normalized.as_str() {
            "full" | "proposed" | "0" => Ok(Self::Full),
            "no-decomposition" | "nodecomposition" | "proposed-nodecomposition" => {
                Ok(Self::NoDecomposition)
            }
            "no-dispatch" | "nodispatch" => Ok(Self::NoDispatch),
            "no-series"
            | "noseries"
            | "no-contraction"
            | "nocontraction"
            | "proposed-nocontraction" => Ok(Self::NoSeries),
            "no-two-cut" | "notwocut" | "no-2-cut" => Ok(Self::NoTwoCut),
            "no-repair" | "norepair" | "proposed-norepair" => Ok(Self::NoRepair),
            "undirected-cubic" | "undirectedcubic" => Ok(Self::UndirectedCubic),
            "one-alternating-worker" | "onealternatingworker" | "one-worker" => {
                Ok(Self::OneAlternatingWorker)
            }
            _ => Err(format!(
                "Unknown ablation '{}'. Expected one of: {}",
                value,
                Self::ALL
                    .iter()
                    .map(|c| c.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aliases_map_to_named_conditions() {
        assert_eq!("Proposed-NoRepair".parse(), Ok(AblationConfig::NoRepair));
        assert_eq!("no_contraction".parse(), Ok(AblationConfig::NoSeries));
        assert_eq!(
            "proposed-nodecomposition".parse(),
            Ok(AblationConfig::NoDecomposition)
        );
    }

    #[test]
    fn ablations_change_only_documented_feature_dependencies() {
        assert!(!AblationConfig::NoDispatch.dispatch_enabled());
        assert!(AblationConfig::NoDispatch.series_enabled());
        assert!(!AblationConfig::NoSeries.alternating_dispatch_enabled());
        assert!(AblationConfig::NoSeries.dispatch_enabled());
        assert!(!AblationConfig::NoDecomposition.two_cut_enabled());
        assert_eq!(
            AblationConfig::OneAlternatingWorker.alternating_workers(),
            1
        );
    }
}
