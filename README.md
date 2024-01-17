# Gtk3NIDAQ

# Installation

1. Clone this repository.
2. cd into `GtkNIDAQ` and clone wsphillip/NIDAQ.jl into the root directory.
3. Launch a repl with the project activated: `julia --project=.`
4. Instantiate the environment: `]instantiate`

When working on the code, always include the project flag (e.g.
`--project=path/to/GtkNIDAQ`) AND specify the threads flag (e.g. `--threads 3,1`; the first
number is number of threads in the default pool, the second number is the number of threads
in the interactive pool.

Example of running a REPL on a POSIX system (use equivalent paths for Windows PC):

```
julia --project=./Gtk3NIDAQ --threads 3,1
# Run code
```
