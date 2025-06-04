import subprocess

def convert_to_cram(
        input_file: str,
        output_file: str,
        threads: int = 1,
        cram_version: str = "3.0",
        write_index: bool = True,
):
    additional_threads = max(0, threads - 1)
    command = [
        "samtools", "view",
        "--output-fmt", f"cram,version={cram_version}",
        "--threads", str(additional_threads),

    ]
    subprocess.run(
    )