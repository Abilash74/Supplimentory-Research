#! /usr/bin/perl -w

($F) = @ARGV;
$J_MAX = 2089;

$i = 0;
open(F, $F) or die "ERROR $F\n";
<F>;
while(<F>){
    chomp;
    @tmp = split(/\t/);
    ($id, $colname, $conc) = ($tmp[0], $tmp[3], $tmp[4]);
    $j = 0;
    if($colname =~/w_(\d+)/){
	$j = $1;
    }else{ next; }
    $id2i{$id} = ++$i if ! $id2i{$id};
    $i2id[$i] = $id if ! $i2id[$i];
    if( $conc == 1){
	$i2label[$i] = "C";
    }else{
	$i2label[$i] = "EO";
    }
    $mat[ $id2i{$id} ][$j] = $conc; 
}
$i_max = $i;

print "ID";
for($j=1; $j<=$J_MAX; $j++){
    print "\tw_" . $j;
}
print "\tlabel\n";
for($i=1; $i<=$i_max; $i++){
    print $i2id[$i];
    for($j=1; $j<=$J_MAX; $j++){
	if( $mat[$i][$j] ){
	    if( $mat[$i][$j] == 1){
		print "\t1";
	    }else{
		printf("\t%.5f", $mat[$i][$j]);
	    }
	}else{
	    print "\t0";
	}
    }
    print "\t" . $i2label[$i] . "\n";
}
