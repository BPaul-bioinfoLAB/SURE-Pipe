#!/usr/bin/env bash
# SURE-Pipe bash completion

_surepipe_complete() {
    local cur prev opts pipeline_params main_nf surepipe_bin surepipe_dir
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # -------------------------------
    # Auto-detect main.nf
    # -------------------------------
    if [[ -n "$SUREPIPE_MAIN_NF" && -f "$SUREPIPE_MAIN_NF" ]]; then
        main_nf="$SUREPIPE_MAIN_NF"
    else
        surepipe_bin="$(command -v SURE-Pipe 2>/dev/null)"
        if [[ -n "$surepipe_bin" ]]; then
            surepipe_dir="$(cd "$(dirname "$surepipe_bin")/.." && pwd)"
            [[ -f "$surepipe_dir/main.nf" ]] && main_nf="$surepipe_dir/main.nf"
        fi

        [[ -z "$main_nf" && -f "./main.nf" ]] && main_nf="./main.nf"
    fi

    [[ -f "$main_nf" ]] || return 0

    # -------------------------------
    # Commands & params
    # -------------------------------
    opts="run --help -h --version -V"

    pipeline_params=$(
        grep -oE 'params\.[A-Za-z0-9_]+' "$main_nf" 2>/dev/null \
        | cut -d. -f2 \
        | sed 's/^/--/' \
        | sort -u
    )

    dir_params="--target_dir --neighbour_dir --output_dir --data_dir --results_dir --reference_dir"
    file_params="--reference_genome --subject_genome --query_genome --input_genomes --genome_pairs"

    # -------------------------------
    # First argument
    # -------------------------------
    if [[ $COMP_CWORD -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "run ${opts}" -- "$cur") )
        return 0
    fi

    # -------------------------------
    # Option name completion
    # -------------------------------
    if [[ "$cur" == --* ]]; then
        COMPREPLY=( $(compgen -W "${pipeline_params} ${opts}" -- "$cur") )
        return 0
    fi

    # -------------------------------
    # Directory completion
    # -------------------------------
    if [[ " $dir_params " == *" $prev "* ]]; then
        COMPREPLY=( $(compgen -d -- "$cur") )
        return 0
    fi

    # -------------------------------
    # File completion
    # -------------------------------
    if [[ " $file_params " == *" $prev "* ]]; then
        COMPREPLY=( $(compgen -f -- "$cur") )
        return 0
    fi

    # -------------------------------
    # --param=value completion
    # -------------------------------
    if [[ "$cur" == --*=* ]]; then
        local pname="${cur%%=*}"
        local pathpart="${cur#*=}"

        if [[ " $dir_params " == *" $pname "* ]]; then
            COMPREPLY=( $(compgen -d -- "$pathpart") )
        else
            COMPREPLY=( $(compgen -f -- "$pathpart") )
        fi

        for i in "${!COMPREPLY[@]}"; do
            COMPREPLY[$i]="${pname}=${COMPREPLY[$i]}"
        done
        return 0
    fi
}

complete -F _surepipe_complete SURE-Pipe
