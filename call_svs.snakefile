configfile: "sv_tools.yaml"
configfile: "data.yaml"

include: "call_svs_sniffles.snakefile"
include: "call_svs_pbsv.snakefile"
include: "call_svs_lumpy.snakefile"
include: "call_svs_grocsvs.snakefile"
include: "call_svs_manta.snakefile"
include: "call_svs_svaba.snakefile"
include: "call_svs_naibr.snakefile"


method_to_sv_call_files_function_dict = {
	"sniffles": expected_sniffles_result_sv_files,
	"pbsv": expected_pbsv_result_sv_files,
	"lumpy": expected_lumpy_result_sv_files,
	"manta": expected_manta_result_sv_files,
	"grocsvs": expected_grocsvs_result_sv_files,
	"svaba": expected_svaba_result_sv_files,
	"naibr": expected_naibr_result_sv_files,
}

def expected_sv_call_files():
	result = []
	for tool_name in config["tools_enabled_methods"]:
		if tool_name in method_to_sv_call_files_function_dict:
			result.extend(method_to_sv_call_files_function_dict[tool_name]())
	return result

rule call_all_svs:
	input: expected_sv_call_files()
