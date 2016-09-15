use ExtUtils::MakeMaker;
my %options = %{       
{
  'TYPEMAPS' => [
    '/usr/share/perl/5.10/ExtUtils/typemap'
  ],
  'NAME' => 'hypergeom_bf98',
  'INC' => '-I/home/pignatelli/workspace/CROC-standAlone',
  'VERSION' => '0.00'
}
};
WriteMakefile(%options);

# Remove the Makefile dependency. Causes problems on a few systems.
sub MY::makefile { '' }
