from ..model import Model
import os, glob, json
from . import comm


class CommitTest(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "committest"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Prepare the test on different abacus commit"
    
    @staticmethod
    def prepare_args(parser):
        parser.add_argument("-j","--jobs",default=[],type=str,help="the path of jobs to be tested",action="extend",nargs="*",)
        parser.add_argument("-i","--image",default="intel",type=str,help="the used image. Can be intel/gnu/cuda that means the latest intel/gnu/cuda image, and also can be a self-defined address. Default is intel", )
        parser.add_argument("-c","--commit",default=[],type=str,help="the commits to be tested. Can also be branch or tag name that can be checkouted to using git checkout",action="extend",nargs="*",)
        parser.add_argument("--repo",default="https://github.com/deepmodeling/abacus-develop.git", help="the repository to be tested. Default is https://github.com/deepmodeling/abacus-develop.git")
        parser.add_argument("--clonemax",default=10,help="the max number of commits to be cloned. Default is 10",type=int)
        parser.add_argument("--install", default="cmake -B build -DENABLE_DEEPKS=ON -DENABLE_LIBXC=ON -DENABLE_LIBRI=ON -DENABLE_RAPIDJSON=ON && cmake --build build -j`nproc` && cmake --install build", help="the command to install the abacus. Default same as dockerfile.intel")
        parser.add_argument("--runcommand",default="OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log", help="the command to run the abacus. Default is OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log")
        parser.add_argument("--machine", default="c32_m64_cpu", help="the machine to run the abacus. Default is c32_m64_cpu")
        parser.add_argument("-r", "--run", default=0, help="if run the test. Default is 0.", type=int)
        
    
    def _scripts(self,maxclone,repo,install,command,script_name):
        scripts = """
pwd=$(pwd)
echo -e "[http]\n       proxy = http://ga.dp.tech:8118\n[https]\n     proxy = http://ga.dp.tech:8118\n" > ~/.gitconfig
"""
        scripts += f"""  
cd /      
max_clone={maxclone}
repo="{repo}"
commit=$1
install="{install}"
command="{command}"
tpath=abacus
clone=0

for i in $(seq 1 $max_clone)
do
    git clone --no-checkout {repo} $tpath
    if [ $? -eq 0 ];then
            clone=1
            break
    fi
    
    echo "Clone failed, retry in 5s"
    sleep 5
done

if [ $clone -eq 0 ];then
    echo "Clone failed"
    exit 1
fi

cd $tpath
git checkout $commit
if [ $? -ne 0 ];then
    echo "Checkout failed"
    exit 1
fi
{install}
if [ $? -ne 0 ];then
    echo "Install failed"
    exit 1  
fi

cd $pwd
{command}
        """
        with open(script_name,"w") as f: f.write(scripts)
        
        return scripts
    
        
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        jobs = params.jobs
        if len(jobs) == 0:
            print("No jobs are specified.")
            return
        commits = params.commit
        if len(commits) == 0:
            print("No commits are specified.")
            return
        
        self._scripts(params.clonemax,params.repo,params.install,params.runcommand,"run_committest.sh")
        
        real_jobs = comm.get_job_list(jobs)
                    
        setting = {
            "save_path": "results",
            "bohrium_group_name": "committest",
            "run_dft": []
        }
        iamge = {"intel": "registry.dp.tech/deepmodeling/abacus-intel:latest",
                 "gnu": "registry.dp.tech/deepmodeling/abacus-gnu:latest",
                 "cuda": "registry.dp.tech/deepmodeling/abacus-cuda:latest"}
        for commit in commits:
            setting["run_dft"].append({
                "sub_save_path": f"{commit}",
                "example": real_jobs,
                "command": f"bash run_committest.sh {commit}",
                "extra_files": ["run_committest.sh"],
                "image": iamge.get(params.image,params.image),
                "bohrium": {
                    "scass_type": params.machine,
                    "job_type": "container",
                    "platform": "ali"
                }
            })
        comm.dump_setting(setting)
        comm.doc_after_prepare("commit_test",real_jobs,["run_committest.sh","setting.json"],has_prepare=False)
        
        #print("After finish the calculation, you can run below command to do the postprocess:")
        
        if params.run:
            bash_script = "ecutwfc_committest.sh"
            with open(bash_script,"w") as f:
                f.write("abacustest submit -p setting.json\n")
                #f.write("cd results\n")
                #f.write(f"abacustest model {self.model_name()} post -m metrics.json\n")
            os.system(f"bash {bash_script} &")
                    
    
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand'''
        pass
    
    
    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        pass