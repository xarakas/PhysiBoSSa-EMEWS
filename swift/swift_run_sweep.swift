import io;
import sys;
import files;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");

app (file out, file err) run_model (file shfile, string param_line, string instance)
{
    "bash" shfile param_line emews_root instance @stdout=out @stderr=err;
}

app (void o) make_dir(string dirname) {
  "mkdir" "-p" dirname;
}

app (void o) cp_config_files(string instance) {
  "cp" "-r" (emews_root+"/data/PhysiBoSSa/config") instance;
}

app (void o) make_output_dir(string instance) {
  "mkdir" "-p" (instance+"/output");
}

file model_sh = input(emews_root+"/scripts/growth_model.sh");
file upf = input(argv("f"));
string upf_lines[] = file_lines(upf);
foreach s,i in upf_lines {
  string instance = "%s/instance_%i/" % (turbine_output, i+1);
  make_dir(instance) => {
    cp_config_files(instance) => {
      make_output_dir(instance) => {
        file out <instance+"out.txt">;
        file err <instance+"err.txt">;
        (out,err) = run_model(model_sh, s, instance);
      }
    }
  }
}