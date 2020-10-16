using Pkg
using Pkg.Artifacts
using Pkg.GitTools

using Pkg.PlatformEngines
Pkg.PlatformEngines.probe_platform_engines!() 

TOML = joinpath(@__DIR__, "..", "notebooks", "Artifacts.toml")
FILES_URL = "https://www.dropbox.com/s/1fd8acem0ivmzcw/files.zip"

files_hash = create_artifact() do artifact_dir
    tarfile = download(FILES_URL)
    @show artifact_dir
    try
        global FILES_HASH = bytes2hex(GitTools.blob_hash(tarfile))
        unpack(tarfile, artifact_dir, verbose=true)
    finally
        rm(tarfile)
    end
end

bind_artifact!(toml, "ESE156_2020", files_hash;
               download_info=[(FILES_URL, FILES_HASH)],
               lazy=true,
               force=true)
