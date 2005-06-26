#!/usr/bin/perl

$pwd = `pwd`;
$pwd =~ s/\n//g;

open PENCILS_INC, "$pwd/src/cparam_pencils.inc";
  @pencils_inc=<PENCILS_INC>;
close PENCILS_INC;

open PENCIL_LIST, "$pwd/pencils.list" or die "Could not find $pwd/pencils.list";
  @pencils=<PENCIL_LIST>;
close PENCIL_LIST;

$npencils=$#pencils;

for ($i=0; $i<$npencils; $i++){
  $pencils[$i] =~ s/\n//g;
  $pencils[$i] =~ s/ //g;
  $pencils[$i] =~ s/%([\dA-Fa-f][\dA-Fa-f])/pack ("C", hex ($1))/eg;
  }

$i=0;

foreach $line (@pencils_inc){
  if ($line =~ /lpencil/) {$pencils_inc[$i]='';}
  $i++
  }

$i=0;
$npencil_def=0;

foreach $line (@pencils_inc){
  if ($line =~ /integer :: i_/){
    $npencil_def++;
    $line_pencil = $line;
    $line_pencil =~ s/^integer :: i_//g;
    $line_pencil =~ s/=.*//g;
    $line_pencil =~ s/\n//g;
    $line_integer = $line;
    $line_integer =~ s/^integer :: i_.*=//g;
    $found=0;
    foreach $pencil (@pencils){
      if ($line_pencil eq $pencil){
        $found=1;
        }
      }
    if ($found == 0) {$lpencil[$i] = '.false.'}
    if ($found == 1) {$lpencil[$i] = '.true.'}
    $i++;
    }
  }

print @pencils_inc;
$lpencil = join ", &\n", @lpencil;
print "logical, parameter, dimension (npencils) :: lpencil=(/ &\n";
print $lpencil;
print "/)"
