
def osg_template(osg_submit,executable,input_files,output_files,arguments,
    cpus= 1,mem= 'GB',Nmem= 1,diskN= 1,diskS= 'GB',log_dir= 'log',queue= 1):

    '''
    Write osg_connect submit file as function.
    - osg_submit: str, submit file path;
    - ID: str, executable path;
    - input_files: str list, all file paths called by program;
    # attention: files must be called by program without path.
    - output_files: str list, output file paths;
    - arguments: str list, flags and values in order.
    '''
    lines= []

    output_dict= {
        x: x.split('/')[-1] for x in output_files
    }

    lines.append('executable = ' + executable)
    lines.append('arguments = ' + ' '.join(arguments))
    lines.append('transfer_input_files = ' +  ','.join(input_files))
    lines.append('transfer_output_files = ' + ','.join(list(output_dict.values())))

    remaps= []
    for filepath,file in output_dict.items():
        remaps.append('='.join([file,filepath]))

    lines.append('transfer_output_remaps = "{}"'.format(' ; '.join(remaps)))

    lines.append('\n')
    lines.append('+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-ubuntu-18.04:latest"')
    lines.append('Requirements = (HAS_MODULES =?= true) && (OSGVO_OS_STRING == "RHEL 7") && (HAS_SINGULARITY == TRUE)')
    lines.append('\n')

    for x in ['error','output','log']:
        lines.append('{} = {}/job.$(Cluster).$(Process).{}.{}'.format(x,log_dir,ID,x))

    lines.append('\n')
    lines.append('request_cpus = ' + str(cpus))
    lines.append('request_memory = {} {}'.format(str(Nmem),mem))
    lines.append('request_disk = {} {}'.format(str(diskN),diskS))

    lines.append('\n')
    lines.append('queue {}'.format(queue))

    with open(osg_submit,'w') as f:
        f.write('\n'.join(lines))
