#!/bin/bash

_surepipe_complete() {
    local cur prev opts pipeline_params main_nf

    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Path of main.nf (auto-detected if user installed properly)
    main_nf="${SUREPIPE_MAIN_NF:-$HOME/SURE-Pipe/main.nf}"

    # CLI commands
    opts="run --help --version -h -V"

    # Auto-extract params.* from Nextflow pipeline
    pipeline_params=$(grep -oP 'params\.\K[A-Za-z0-9_]+' "$main_nf" | sed 's/^/--/')

    # Autocomplete parameters beginning with --
    if [[ ${cur} == --* ]]; then
        COMPREPLY=( $(compgen -W "${pipeline_params}" -- "${cur}") )
        return 0
    fi

    # After 'SURE-Pipe'
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
        return 0
    fi

    # After 'run'
    if [[ "${prev}" == "run" ]]; then
        COMPREPLY=( $(compgen -W "${pipeline_params}" -- "${cur}") )
        return 0
    fi
}

complete -F _surepipe_complete SURE-Pipe
