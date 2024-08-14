project('mpiabc', 'c')

gsl_dep = dependency('gsl')
libabc = library('abc', 'model.c', 'sampler.c',
                 dependencies : gsl_dep)
executable('infer', 'infer.c',
           link_with : libabc,
           dependencies : gsl_dep)

check_dep = dependency('check',
                       fallback : ['check', 'check_dep'],
                       # Suppress -Wformat-extra-args messages.
                       default_options : ['warning_level=0'])
text_exe = executable('test_abc', 'test_abc.c',
                      link_with : libabc,
                      dependencies : [gsl_dep, check_dep])
test('./test_abc', text_exe)

conf_data = configuration_data()
conf_data.set('INCLUDE_PWD', meson.project_source_root())
conf_data.set('INCLUDE_GSL', gsl_dep.get_variable('includedir'))
if check_dep.type_name() == 'pkgconfig'
  conf_data.set('INCLUDE_CHECK', check_dep.get_variable('includedir'))
else
  conf_data.set('INCLUDE_CHECK',
                meson.current_build_dir() / 'subprojects' /
                'check' + '-' + check_dep.version() / 'src')
endif
configure_file(input : '.dir-locals.el.in',
               output : '.dir-locals.el',
               configuration : conf_data)
