let knm = Array.init 3 (function i -> Array.create 3 0.0)
;; 
let fill_array () =
    Printf.fprintf stderr "filling\n%!";
    for k = 0 to 2 do
    for f = 0 to 2 do
        Printf.fprintf stderr "%d %d\n%!" k f;
        let p = (float_of_int f /. 100.0) in
        knm.(k).(f) <-  5.0 ; 
    done
 done;
 Printf.fprintf stderr "end of filling\n"
;;
let _ = ignore (fill_array ())
;;

