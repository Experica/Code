## Lossless encode each epoch of images to a video
using FFMPEG

function imagevideogroup(dir;iext=".TIFF",vext=".mkv")
    fns = filter!(f->endswith(f,iext),readdir(dir))
    ns = map(f->match(Regex("(\\w*)-\\d*$iext"),f).captures[1],fns)
    map(n->(;images=joinpath(dir,"$n-%d$iext"),video=joinpath(dir,"$n$vext")),unique!(ns))
end



for g in imagevideogroup("C:\\Users\\fff00\\Pictures\\Test")
    @ffmpeg_cmd `-y -i $(g.images) -pix_fmt:v + -c:v ffv1 -level 3 $(g.video)`
end



@ffmpeg_cmd `-y -i $images -pix_fmt:v gray16le -c:v jpeg2000 -format j2k -pred 1 $video.avi`
@ffmpeg_cmd `-y -i $images -pix_fmt:v gray16le -c:v libx264 -preset ultrafast -crf 0 $video.mp4`
@ffmpeg_cmd `-y -i $images -pix_fmt:v gray16le -c:v libx265 -preset ultrafast -x265-params lossless=1 $video.mp4`

FFMPEG.exe(`-pix_fmts`)

## Test frame lossless
using VideoIO,Images

fs = VideoIO.load("C:\\Users\\fff00\\Pictures\\Test\\Test.mkv")



vh = openvideo("C:\\Users\\fff00\\Pictures\\Test1\\Test1.mkv")

v0=read(vh)

t=v[1][1,1]

N0f16(t)
reinterpret(UInt16,collect(v0))
close(vh)

f1 = load("C:\\Users\\fff00\\Pictures\\Test\\Test-0.TIFF")



VideoIO.VIO_PIX_FMT_DEF_ELTYPE_LU
