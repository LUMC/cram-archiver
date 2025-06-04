import subprocess

# 3.1 Not supported by some tools currently (2025)
DEFAULT_CRAM_VERSION = "3.0"



def convert_to_cram(
        input_file: str,
        output_file: str,
        reference: str,
        threads: int = 1,
        cram_version: str = DEFAULT_CRAM_VERSION,
        write_index: bool = True,
):
    additional_threads = max(0, threads - 1)
    command = [
        "samtools", "view",
        "--output-fmt", f"cram,version={cram_version}",
        "--threads", str(additional_threads),
        "-o", output_file,
        "--reference", reference,
        input_file,
    ]
    if write_index:
        command.append("--write-index")
    subprocess.run(command, check=True)
