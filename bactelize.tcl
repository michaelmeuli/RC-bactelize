#!/usr/bin/tclsh


package require json
set infile [open "output.json" r]
set file_data [read $infile]
puts $file_data
set parsed [json::json2dict $file_data]
puts "\n\n"




dict for {k1 v1} $parsed {
    puts "$k1 => $v1"
    puts "\n\n"
    # Iterate over contents
    dict for {k2 v2} $v1 {
        puts "    $k2 => $v2"
        puts "\n\n"
        dict for {k3 v3} $v2 {
            puts "        $k3 => $v3"
            puts "\n\n"
	}

    }
    # Picking out a particular key
    #puts "The grade was [dict get $val grade]" 
}
puts "\n\n"
puts [dict get $parsed params init init_blob_min]
dict set parsed params init init_blob_min 20
puts [dict get $parsed params init init_blob_min]





package require Tcl 8.6
package require json::write
 
proc tcl2json value {
    # Guess the type of the value; deep *UNSUPPORTED* magic!
    regexp {^value is a (.*?) with a refcount} \
	[::tcl::unsupported::representation $value] -> type
    
    puts $type

    switch $type {
	string {
	    return [json::write string $value]
	}
	dict {
	    return [json::write object {*}[
		dict map {k v} $value {tcl2json $v}]]
	}
	list {
	    return [json::write array {*}[lmap v $value {tcl2json $v}]]
	}
	int - double {
	    return [expr {$value}]
	}
	booleanString {
	    return [expr {$value ? "true" : "false"}]
	}
	default {
	    # Some other type; do some guessing...
	    if {$value eq "null"} {
		# Tcl has *no* null value at all; empty strings are semantically
		# different and absent variables aren't values. So cheat!
		return $value
	    } elseif {[string is integer -strict $value]} {
		return [expr {$value}]
	    } elseif {[string is double -strict $value]} {
		return [expr {$value}]
	    } elseif {[string is boolean -strict $value]} {
		return [expr {$value ? "true" : "false"}]
	    }
	    return [json::write string $value]
	}
    }
}

#set d [dict create blue [list 1 2] ocean water]
#puts [tcl2json $parsed]

set fp [open "output2.json" w+]
puts $fp [tcl2json $parsed]
close $fp


