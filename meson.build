project('mpiabc', 'c')

# Target: infer
gsl_dep = dependency('gsl')
libabc = library('abc', 'model.c', 'sampler.c',
                 dependencies : gsl_dep)
executable('infer', 'infer.c',
           link_with : libabc,
           dependencies : gsl_dep)

# Optional target: latex/refman.pdf
doxygen = find_program('doxygen', version : '>=1.9.8', required : false)
sed = find_program('sed', required : false)
gmake = find_program('gmake', required : false)
if doxygen.found() and sed.found() and gmake.found()
  message('Optional target "docs" possible that uses doxygen, sed, and gmake')

  doxyfile_data = configuration_data()
  doxyfile_data.set('INPUT', meson.project_source_root())
  doxyfile = configure_file(input : 'Doxyfile.in',
                            output : 'Doxyfile',
                            configuration : doxyfile_data,
                            install : false)
  docs_latex = custom_target('docs-latex',
                             command: [meson.project_source_root() /
                                       'doxygen.sh',
                                       doxyfile],
                             input : [
                               'abc.h',
                               'citations.bib',
                               'infer.c',
                               'model.c',
                               'model.h',
                               'README.md',
                               'sampler.c',
                               'sampler.h',
                             ],
                             output : 'latex',
                             build_by_default : false)
  run_target('docs',
             command : [gmake, '-C', 'latex'],
             depends : [docs_latex])
else
  warning('Optional target "docs" not possible because it requires doxygen, sed, and gmake')
endif

# Optional target: test_abc
check_dep = dependency('check',
                       fallback : ['check', 'check_dep'],
                       # Suppress -Wformat-extra-args messages.
                       default_options : ['warning_level=0'])
text_exe = executable('test_abc', 'test_abc.c',
                      link_with : libabc,
                      dependencies : [gsl_dep, check_dep])
test('test_abc', text_exe)

# Target: .dir_locals.el
dir_locals_data = configuration_data()
dir_locals_data.set('INCLUDE_PWD', meson.project_source_root())
dir_locals_data.set('INCLUDE_GSL', gsl_dep.get_variable('includedir'))
if check_dep.type_name() == 'pkgconfig'
  dir_locals_data.set('INCLUDE_CHECK', check_dep.get_variable('includedir'))
else
  dir_locals_data.set('INCLUDE_CHECK',
                      meson.current_build_dir() / 'subprojects' /
                      'check' + '-' + check_dep.version() / 'src')
endif
configure_file(input : '.dir-locals.el.in',
               output : '.dir-locals.el',
               configuration : dir_locals_data)
