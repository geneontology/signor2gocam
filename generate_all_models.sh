mkdir -p $2
for f in $1/* ; do
  f_base=$(basename $f)
  f_base=${f_base%.tsv}
  title=$(cut -f2 $f | head -n 2 | tail -n 1)
  cmd="python3 pathway_importer.py -f $f -t SIGNOR - $title -o $2/$f_base.ttl"
  echo $cmd
  $cmd
done