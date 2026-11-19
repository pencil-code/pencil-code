import os
import re
import sys
import argparse

CONFIG_FILES = "GNU-GCC_MPI+GNU-GCC_GPU+GNU-GCC_debug"


def gpu_autotest_samples(pc_home):
    """Return the samples listed in @gpu_tests of bin/pc_auto-test (commented-out entries skipped)."""
    with open(os.path.join(pc_home, "bin", "pc_auto-test")) as f:
        text = f.read()
    match = re.search(r"my\s+@gpu_tests\s*=\s*qw\((.*?)^\s*\);", text, re.S | re.M)
    if not match:
        raise RuntimeError("could not find @gpu_tests in bin/pc_auto-test")
    samples = []
    for line in match.group(1).splitlines():
        for sample in line.split("#", 1)[0].split():
            if sample not in samples:
                samples.append(sample)
    return samples


def sample_dir(pc_home, sample):
    """Resolve a sample name: absolute path, path relative to samples/, or a build-samples/<name>."""
    if os.path.isabs(sample):
        return sample
    for candidate in (os.path.join(pc_home, "samples", sample),
                      os.path.join(pc_home, "samples", "build-samples", sample)):
        if os.path.isdir(candidate):
            return candidate
    raise RuntimeError(f"sample not found: {sample}")


def run(command):
    print(f"+ {command}", flush=True)
    return os.system(command)


def update_sample(pc_home, sample, verify):
    directory = sample_dir(pc_home, sample)
    print(f"=== {sample} ({directory})", flush=True)
    build = (f"cd {directory} && pc_setupsrc --force-astaroth && "
             f"pc_build -f {CONFIG_FILES} MODIFY_SOURCE_CODE=on")
    status = run(build)
    if status or not verify:
        return status
    # check that the modified sources still build from scratch
    rebuild = (f"cd {directory} && make distclean ; rm -rf src/astaroth && "
               f"pc_setupsrc --force-astaroth && pc_build -f {CONFIG_FILES} FFLAGS+=-O0")
    return run(rebuild)


def main():
    argparser = argparse.ArgumentParser(description="Update *.f90 files to be GPU ready",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("-p", "--pencil-home", help="Where PC is installed", required=True)
    group = argparser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--sample", action="append",
                       help="Sample to run over (relative to samples/, or a build-samples name); may be repeated")
    group.add_argument("-g", "--gpu-autotest", action="store_true",
                       help="Run over all samples of the GPU auto-test (@gpu_tests in bin/pc_auto-test)")
    argparser.add_argument("-l", "--list", action="store_true", help="Only print the samples and exit")
    argparser.add_argument("-v", "--verify", action="store_true",
                           help="After updating a sample, rebuild it from scratch to check the modified code compiles")
    args = argparser.parse_args()
    pc_home = args.pencil_home
    os.environ["PENCIL_USER"] = "GPU_AUTOTEST"

    samples = gpu_autotest_samples(pc_home) if args.gpu_autotest else args.sample
    if args.list:
        print("\n".join(samples))
        return 0

    status = run("git submodule update --init --remote")
    if status:
        return 1
    for sample in samples:
        # stop at the first failure so that broken code does not get committed
        if update_sample(pc_home, sample, args.verify):
            print(f"FAILED: {sample}", flush=True)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
