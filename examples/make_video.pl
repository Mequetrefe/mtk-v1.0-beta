#!/usr/bin/perl

$nframes = 9;
$delay = 10;
$name = "anim-ex09";

for ($j=1; $j<=$nframes; $j++) {
  $list .= "snsp_"."$j"."00.png ";
}
`convert -delay $delay $list $name.gif`;
`rm -f *.dat *.gpl snsp*`
