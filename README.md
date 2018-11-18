## Bio Algorithms

This repository contains various algorithms and data structures used to solve the assignments for the course "Algorithmen der Bioinformatik" (Algorithms of bioinformatics) at the University of Göttingen.

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
After the Rust toolchain has been installed, the project can be build by entering
```
cargo build --release
```
in the project root.  
The documentation can be built and opened by entering
```
cargo doc --no-deps --open
```
