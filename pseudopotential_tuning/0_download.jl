using ASEconvert

download_url = [
    "https://raw.githubusercontent.com/aiidateam/acwf-verification-scripts/main/0-preliminary-do-not-run/unaries/xsfs-unaries-verification-PBE-v1/Li-BCC.xsf",
    "https://raw.githubusercontent.com/aiidateam/acwf-verification-scripts/main/0-preliminary-do-not-run/oxides/xsfs-oxides-verification-PBE-v1/Li-XO.xsf",
]


tmp_dir = tempdir()
out_dir = @__DIR__

for url in download_url
    @info "Downloading:" url
    bname = basename(url)
    file = download(url, joinpath(tmp_dir, bname))
    
    # Convert from xsf to extxyz
    name, ext = splitext(bname)
    system = ase.io.read(file)
    system.info["name"] = name
    ase.io.write(joinpath(out_dir, "$name.extxyz"), system)
end

