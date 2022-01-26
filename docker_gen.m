appFile = "SbS_indexesFromEMGinEBformat.m"
%buildResults = compiler.build.standaloneApplication(appFile, 'AdditionalFiles', "tools/read_simple_yaml.m")
buildResults = compiler.build.standaloneApplication(appFile)

opts = compiler.package.DockerOptions(buildResults, 'ImageName','pi_sbs_emg')

compiler.package.docker(buildResults, 'Options', opts)
