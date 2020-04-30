mkdir -p $2
for f in $1/* ; do
  f_base=$(basename $f)
  f_base=${f_base%.tsv}
  cmd="python3 pathway_importer.py -f $f -t $f_base -o $2/$f_base.ttl"
  echo $cmd
  $cmd
done