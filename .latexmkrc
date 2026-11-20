# Run latexmk from the repository root.  Main documents live directly in papers/;
# put fragments in papers/sections/.  No filename edit is needed after renaming
# paper-template.tex to <repository-name>.tex.
@default_files = ('papers/*.tex');

# Resolve \input and bibliography paths from each main document's directory.
$do_cd = 1;
$pdf_mode = 4;  # LuaLaTeX
$lualatex = 'lualatex -synctex=1 -interaction=nonstopmode -halt-on-error -file-line-error %O %S';
$bibtex_use = 2;

# With main documents in papers/ and $do_cd=1, this is the repository's build/.
$out_dir = '../build';

# Keep SyncTeX beside the PDF for the editor's forward and inverse search.
# Do not enable continuous preview (-pvc); VS Code can build when saving.
