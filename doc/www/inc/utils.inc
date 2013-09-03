<?
	# $Id$

        # Echoes "$pre$str$post$always" or just "$always$default" if $str is empty.
        function ifecho ($pre, $str, $post, $always = "", $default = "") {
                $str = ifstr ($pre, $str, $post, $always, $default);
                if (!empty ($str)) { echo $str; }
        }

        # Returns "$pre$str$post$always" or just "$always$default" if $str is empty.
        function ifstr ($pre, $str, $post, $always = "", $default = "") {
                if (isset ($str) && is_string ($str)) {
                        $str = trim ($str);
                        if (!empty ($str)) {
                                return ($pre.$str.$post.$always);
                        }
                }
                return ($always.$default);
        }

        # Echoes a link to $link labeled $str if $link is not empty.
        function iflink ($link, $str, $class = "", $target = "", $name = "", $js = "") {
                if ($link) { echolink ($link, $str, $class, $target, $name, $js); }
        }

        # Returns a link to $link labeled $str or just $str if $link is empty.
        function genlink ($link, $str, $class = "", $target = "", $name = "", $js = "") {
                $link = str_replace (array ("&", " "), array ("&amp;", "+"), $link);
                $attributes = ifstr (' NAME="', $name, '"').ifstr (' CLASS="', $class, '"').ifstr (' TARGET="', $target, '"').ifstr (" ", $js, "");
                return (ifstr ('<A HREF="', $link, '"'.$attributes.'>'.$str.'</A>', "", $str));
        }

        # Echoes a link to $link labeled $str or just $str if $link is empty.
        function echolink ($link, $str, $class = "", $target = "", $name = "", $js = "") {
                echo genlink ($link, $str, $class, $target, $name, $js);
        }
 ?>

