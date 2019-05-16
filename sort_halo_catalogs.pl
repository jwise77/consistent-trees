#!/usr/bin/perl -w
use lib qw(src);
use ReadConfig;
use MultiThread;
use Universe::Time;
use MassConversions;
use IO::File;
use Scalar::Util;

my $config = new ReadConfig($ARGV[0]);
    $config->set_defaults(
        HLIST_OUTBASE => "/Volumes/Peter 2/Bolshoi/Hlists");
our $HLIST_OUTBASE = $config->{HLIST_OUTBASE};


$MultiThread::threadlimit = 1;

opendir DIR, $HLIST_OUTBASE;
my @hlists = grep { /^hlist.*\.list$/ } readdir DIR;
closedir DIR;

for (@hlists) {
    my $fn = "$HLIST_OUTBASE/$_";
    next unless MultiThread::ForkChild();
    sort_hlist($fn);
    exit;
}
MultiThread::WaitForChildren();

sub sort_hlist {
    my $fn = shift;
    system("./sort_halo_catalogs", $fn);
}


sub sort_halos {
    my (undef, $upid_a, $m_a, $hm_a) = @$a;
    my (undef, $upid_b, $m_b, $hm_b) = @$b;
    return -1 if ($hm_a > $hm_b);
    return 1 if ($hm_b > $hm_a);
    return -1 if ($upid_a < $upid_b);
    return 1 if ($upid_a > $upid_b);
    return ($m_b <=> $m_a);
}

sub sort_hlist_old {
    my $fn = shift;
    my @data;
    open INPUT, "<", $fn;
    open OUTPUT, ">", "$fn.sorted";
    while (<INPUT>) {
	if (/^#/) { print OUTPUT; next; }
	my @a = split;
	my $i = @data;
	$hosts{$a[1]} = $a[61];
	if ($a[6] < 0) { $a[6] = $a[1]; }
	push @data, [$i, $a[6], $a[61]];
	push @halos, $_;
	$upids{$a[1]} = $a[6];
    }
    for (@data) {
	while ($upids{$_->[1]} > -1) { $_->[1] = $upids{$_->[1]}; }
	$_->[3] = $hosts{$_->[1]};
	print "@$_\n" unless ($hosts{$_->[1]});
    }
    @data = sort sort_halos @data;
    for (@data) {
	print OUTPUT $halos[$_->[0]];
    }
    close OUTPUT;
}


