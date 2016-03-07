#!/usr/bin/perl -w

# prog1.1 
# Bruno Contreras-Moreira
# Nearest Neighbor dG calculator

#use strict;
use Data::Dumper;

# global variables
my $T           = 37; # temperature(C)
my $windowL_E1  = 50;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my $windowL_E2	= 100; #window length 2
my $windowL	= 15; #tamano de la ventana en nucleotidos
my $window_dist	= 50; #distance between windows	
my %Secuencias_K12;
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/k�mol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H',   0, 'S',-1.4 } );

my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";

print "Parameters: Temperature=$T\\C Window=$windowL_E1, $window_dist, $windowL_E2\n\n";

open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>)
{
	##Almacenar el nombre y la secuencia del archivo
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		my ($name,$seq) = ($1,$2); 
		my $len_seq=length($seq);
		print("sequence $name : $len_seq nts \n");
		$Secuencias_K12{$name}=$seq;
	}
}
close(SEQ);


# calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
# parameters: 1) DNA sequence string; 2) Celsius temperature
# returns; 1)  Hash con tres llaves: 1)Vector que contiene el valor de D(n) para cada n de la secuencia que esté dentro del rango permitido 2)Valores de E1(n), 3)Valores de E2(n).
# uses global hash %NNparams
# uses global hash %Secuencias_K12
sub duplex_deltaG
{
   	my ($seq,$tCelsius) = @_; 

	my @dG_individuales;
	my @free_energies;
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);
	my @sequence = split(//,uc($seq));
	my $tK = 273.15 + $tCelsius;
	
	sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }
	
	# add dG for overlapping dinculeotides
	for(my $i=0;$i<$#sequence;$i++) 
	{			
			$DNAstep = $sequence[$i].$sequence[$i+1].'/'.
				complement($sequence[$i].$sequence[$i+1]);
			
			if(!defined($NNparams{$DNAstep}))
			{
				$DNAstep = reverse($DNAstep);
			}
			
			$dG = ((1000*$NNparams{$DNAstep}{'H'})-
					($tK*$NNparams{$DNAstep}{'S'}))
					/ 1000;
			push @dG_individuales, $dG; #dg_individuales es un vectr donde se almacena la deltaG entre cada par de nucleotidos
	}
	
	##Calcular el numero de ventaas y de n's posibles
	my $numero_vent=(scalar @sequence)-$windowL+1;
	my $numero_n= $numero_vent-($windowL_E1+$windowL_E2+$window_dist-1);
	my @E1;
	my @E2;
	my @D;

	#Sacar la energia libre para cada n con una ventana de 15 nt.
	for(my $i=0; $i<$numero_vent; $i++){
		$free_energies[$i]=0;
		for(my $j=0; $j<$windowL; $j++){
			$free_energies[$i]+=$dG_individuales[$i+$j];
		}
	}

	# add correction for helix initiation
	for(my $i=0; $i<$numero_vent; $i++){	
		$nt = $sequence[$i]; # first nucleotide of the window
		if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
		$free_energies[$i] += ((1000*$NNparams{$nt}{'H'})-
						($tK*$NNparams{$nt}{'S'}))
						/ 1000; 
		$nt = $sequence[$i+$windowL-1]; # last nucleotide of the window
		if(!defined($NNparams{$nt})){ $nt = complement($nt) }
		$free_energies[$i] += ((1000*$NNparams{$nt}{'H'})-
						($tK*$NNparams{$nt}{'S'}))
						/ 1000;
	}				
	# symmetry correction
	for(my $i=0; $i<$numero_vent; $i++){
		my $seq_ventana;
		for(my $j=0; $j<$windowL; $j++){
			$seq_ventana=$seq_ventana.$sequence[$i+$j];
		}
		if($seq_ventana eq reverse($seq_ventana)){
			$free_energies[$i] += ((1000*$NNparams{'sym'}{'H'})-
					($tK*$NNparams{'sym'}{'S'}))
					/ 1000;
		}
	}

   #Sacar el valor de E1 para cada n
	for(my $i=0; $i<$numero_n; $i++){
		$E1[$i]=0;
		for(my $j=0; $j<$windowL_E1; $j++){
			$E1[$i]+=$free_energies[$i+$j];
		}
	}
	#Calculate the average
	for(my $i=0; $i<$numero_n; $i++){
		$E1[$i]=$E1[$i]/$windowL_E1;
	}
	
    #Sacar el valor de E2 para cada n
	for(my $i=0; $i<$numero_n; $i++){
		$E2[$i]=0;
		for(my $j=0; $j<$windowL_E2; $j++){
			$E2[$i]+=$free_energies[$windowL_E1+$window_dist+$i+$j];
		}
	}

	#Calculate the average
	for(my $i=0; $i<$numero_n; $i++){
		$E2[$i]=$E2[$i]/$windowL_E2;
	}

   #Sacar D(n)
	for (my $i=0; $i<$numero_n; $i++){
		$D[$i]=$E1[$i]-$E2[$i];
	}

	my %D_E;
	$D_E{'D'}=\@D;
	$D_E{'E1'}=\@E1;
	$D_E{'E2'}=\@E2;
	return %D_E;
}




my $llave;
my %deltaGs;
my @D_n;
my @E1_n;
my @E2_n;

my $filename = 'Prediccion_promotores.txt';

open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
foreach $llave (sort(keys %Secuencias_K12)){
	my %result=duplex_deltaG($Secuencias_K12{$llave}, $T);

	$deltaGs{$llave}= \%result;

	my $seq=$Secuencias_K12{$llave};
	my $numero_n=(length $seq)-($windowL+$windowL_E1+$windowL_E2+$window_dist-2);
	print $fh ("Nombre: $llave\n");

	for (my $i=0; $i<$numero_n; $i++){
		$D_n[$i]=0;
		$E1_n[$i]=0;
		$E2_n[$i]=0;
	}

	for(my $i=0; $i<$numero_n; $i++){
		print $fh ("\tD(n) $i: $result{'D'}[$i]\n");
		$D_n[$i]+=$result{'D'}[$i];
		
	}
	for(my $i=0; $i<$numero_n; $i++){
		print $fh ("\tE1(n) $i: $result{'E1'}[$i]\n");
		$E1_n[$i]+=$result{'E1'}[$i];
		$E2_n[$i]+=$result{'E2'}[$i];
	}
}
close $fh;
print "done\n";


###########################

#Calcular el promedio de D, E1 y E2 para cada n

my $filename = 'Promedio_promotores.txt';
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print $fh "D(n)\n";
print $fh join(", ", @D_n);
print $fh "\nE1(n)\n";
print $fh join(", ", @E1_n);
print $fh "\nE2(n)\n";
print $fh join(", ", @E2_n);
close $fh;
print "done\n";

