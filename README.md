## Bio Algorithms

This repository contains various algorithms and data structures used to solve the assignments for the course "Algorithmen der Bioinformatik" (Algorithms of bioinformatics) at the University of Göttingen.  
The current docs can be found [here](https://robinwilliam.hundt.pages.gwdg.de/bio_algos/bio_algos/index.html).  
The latest binary compiled for Linux can be downloaded from [here](https://gitlab.gwdg.de/robinwilliam.hundt/bio_algos/-/jobs/artifacts/master/file/target/release/bio_algos?job=build_release). 

### Why Rust?
The [Rust](https://www.rust-lang.org/en-US/index.html) programming language is "a systems programming language that runs blazingly fast, prevents segfaults, and guarantees thread safety."  
It features:
- zero-cost abstractions
- move semantics
- guaranteed memory safety
- threads without data races
- trait-based generics
- pattern matching
- type inference
- minimal runtime
- efficient C bindings

This makes it an ideal match for implementing efficient and correct algorithms and data structures.  
The strong safety guarantees made by the language allow the easy parallelization of even complex algorithms (at the moment parallelization is not used in this crate). 


### Installing Rust and building the project
The recommended way of installing Rust is via the `rustup` tool. See the [Install Rust](https://www.rust-lang.org/en-US/install.html) page for installation instructions.  
This project requires features that are only available in the nightly toolchain of Rust, it can be installed via rustup:
```
rustup install nightly
rustup default nightly
```
After the Rust nightly toolchain has been installed, the project can be build by entering
```
cargo build --release
```
in the project root. The binary will be emitted at `target/release/bio_algos`.  
The documentation can be built and opened by entering
```
cargo doc --no-deps --open
```

The few tests specified in the analysis module can be run by entering:
```
cargo test
```

### Protocol
The current version of the assignment can, including compiled binary and execute + to html converted notebook can be found [here](https://gitlab.gwdg.de/robinwilliam.hundt/bio_algos/-/jobs/artifacts/master/file/hundt_robin.tar.gz?job=build_assignment).  
The html version of the protocol will be under `pwm-assignment/protocol.html`.