use log::LevelFilter;

use tosa::cli;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse CLI arguments
    let matches = cli::build_cli().get_matches();
    let config = cli::parse_config(&matches);

    // Initialize logger
    if config.verbose {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Info)
            .init();
    }

    // Run the pipeline
    tosa::run(&config)
}
