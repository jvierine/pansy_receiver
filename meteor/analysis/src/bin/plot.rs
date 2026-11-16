use anyhow::{ensure, Result};
use std::path::Path;

fn main() -> Result<()> {
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() == 4 && args[1] == "radiants" {
        return meteor_analysis::plotting::radiant_figure(Path::new(&args[2]), Path::new(&args[3]));
    }
    ensure!(
        args.len() == 3,
        "Usage: meteor-plot ANALYSIS_DIRECTORY FIGURE_DIRECTORY"
    );
    meteor_analysis::plotting::plot_run(Path::new(&args[1]), Path::new(&args[2]))
}
