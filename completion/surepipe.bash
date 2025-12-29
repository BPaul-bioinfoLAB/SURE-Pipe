#!/usr/bin/env bash
# SURE-Pipe bash completion (final, ls-like)

_surepipe_complete() {
    local cur prev
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # -------------------------------
    # Static options
    # -------------------------------
    local opts="run --help -h --version -V"

    local dir_params="--target_dir --neighbour_dir --output_dir --data_dir --results_dir --reference_dir"
    local file_params="--reference_genome --subject_genome --query_genome --input_genomes --genome_pairs --target_genome"

    # -------------------------------
    # First argument
    # -------------------------------
    if [[ $COMP_CWORD -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "run ${opts}" -- "$cur") )
        return 0
    fi

    # -------------------------------
    # Directory params → directory completion
    # -------------------------------
    if [[ " $dir_params " == *" $prev "* ]]; then
        _filedir -d
        return 0
    fi

    # -------------------------------
    # File params → ls-style file completion
    # THIS IS THE KEY FIX
    # -------------------------------
    if [[ " $file_params " == *" $prev "* ]]; then
        _filedir
        return 0
    fi

    # -------------------------------
    # --param=value form
    # -------------------------------
    if [[ "$cur" == --*=* ]]; then
        local pname="${cur%%=*}"

        if [[ " $dir_params " == *" $pname "* ]]; then
            _filedir -d
        else
            _filedir
        fi

        # Reattach --param=
        for i in "${!COMPREPLY[@]}"; do
            COMPREPLY[$i]="${pname}=${COMPREPLY[$i]}"
        done
        return 0
    fi

    # -------------------------------
    # Option name completion
    # -------------------------------
    if [[ "$cur" == --* ]]; then
        COMPREPLY=( $(compgen -W "${opts}" -- "$cur") )
        return 0
    fi
}

complete -F _surepipe_complete SURE-Pipe
