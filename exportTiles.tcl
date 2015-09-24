# Function to extract all tiles from a LIF file
# The files are all saved as separate files in the destination directory.
#
# Note: the destination directory is optional. If it is not given the current
# working directory will be used.
proc extractTilesFromLif {filename {dir {.}}} {
    set metadata    [img preOpen $filename]
    set mainCounter 0

    dict for {key value} [dict get $metadata subImages] {
        set subCounter 0

        foreach tile [dict get $value subImages] {
            huOpt report "Reading: $key - $tile"
            
            set img [img open $filename \
                         -subImage $key \
                         -subIndex $mainCounter:$subCounter \
                         -tileIndex $subCounter]
            $img show
            $img save [file join $dir "export_${key}_${tile}"]
            $img del
            
            incr subCounter
        }
        incr mainCounter
    }
}
