#  bash autocompletion
if [ -n "$BASH" ]
then
	 _nvg_complete()
	 {
		  local cur prev
		  COMPREPLY=()
		  cur="${COMP_WORDS[COMP_CWORD]}"
		 opts="align annotate augment autoindex call chains chunk circularize clip combine construct convert deconstruct depth describe dotplot filter find gamcompare gampcompare gamsort gbwt giraffe haplotypes help ids index inject kmers map mask minimizer mod mpmap pack paths primers prune rna sim simplify snarls stats surject trace validate vectorize version view"

		  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
		  return 0
	 }
	 complete -o nospace -F _nvg_complete nvg
fi
