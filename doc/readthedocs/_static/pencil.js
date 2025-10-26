var Scorer = {
    // Implement the following function to further tweak the score for
    // each result The function takes a result array [filename, title,
    // anchor, descr, score] and returns the new score.
    score: function (result) {
        let updated = result[4];
        // Push code documentation down so "regular" docs will show
        // up before them.
        if (result[0].search("code/cpp/") >= 0) {
            updated -= 15;
        } else if (result[0].search("code/py/") >= 0) {
            updated -= 14;
        } else if (result[0].search("code/js/") >= 0) {
            updated -= 13;
        } else if (result[0] == "index") {
            // Index file lists everything; that's probably not very helpful
            updated -= 3;
        }

        // The deeper the document is nested, the more it should be penalized
        if (result[0].search("code/") >= 0) {
            updated -= 2 * (result[0].split("/").length - 1);
        }
        return updated;
    },

    // query matches the full name of an object
    objNameMatch: 7,
    // or matches in the last dotted part of the object name
    objPartialMatch: 6,
    // Additive scores depending on the priority of the object
    objPrio: {
        0: 15,  // used to be importantResults
        1: 5,   // used to be objectResults
        2: -5
    },  // used to be unimportantResults
    //  Used when the priority is not in the mapping.
    objPrioDefault: 0,

    // query found in title
    title: 8,
    partialTitle: 5,
    // query found in terms
    term: 2,
    partialTerm: 1
};

$(document).ready(function () {
    $(".fa-bitbucket").each(function (e) {
        /* Display only after changing text */
        if ($(this).html() === " Edit on Bitbucket") {
            $(this).html(" View on Bitbucket");
        }
        $(this).addClass("fa-bitbucket-show");
    });
    $(".collapsable ul").each(function (i) {
        if ($(this).children().length <= 5) return;
        $(this).parent().addClass("closed")
        $(this).parent().append('<div class="control"></div>')
        let ctrl = $(this).siblings(".control");
        ctrl.html('<a href="#">Show more <i class="fa fa-plus"></i></a>');
        ctrl.css("margin-top", "5px");
        let margin = $(this).css("margin-bottom");
        $(this).css("margin-bottom", "0");
        $(this).parent().css("margin-bottom", margin);
    });
    $(".control").click(function (e) {
        let parent = $(this).parent();
        if (parent.is(".closed")) {
            $(this).html('<a href="#">Show less <i class="fa fa-minus"></i></a>');
            var isClosing = false;
        } else {
            $(this).html('<a href="#">Show more <i class="fa fa-plus"></i></a>');
            var isClosing = true;
        }
        $(this).parent().toggleClass("closed");
        if (isClosing) {
            $(this).parent().prev("p").scrollintoview();
        }
        return false;
    });
    $("#customvarsearch").keyup($.debounce(250, function () {
        const needle = $(this).val().toLowerCase();
        let hidden = 0, shown = 0;
        if (!needle) {
            $('section[id^=module-] table tbody').children("tr").each(function () {
                $(this).show();
            });
            $('section[id^=module-]').each(function (idx, el) {
                $(this).show();
            });
        } else {
            $('section[id^=module-] table tbody').children("tr").each(function () {
                let found = false;
                $(this).children("td:first-child").each(function (idx, el) {
                    if ($(this).text().toLowerCase().indexOf(needle) > -1) {
                        found = true;
                    }
                });
                if (found) {
                    $(this).show();
                    $(this).parents('section[id^=module-]').show();
                    shown += 1;
                } else {
                    $(this).hide();
                    hidden += 1;
                }
            });
            $('section[id^=module-]').each(function (idx, el) {
                if ($(`#${el.id} table tbody tr td :visible`).length === 0) {
                    $(this).hide();
                } else {
                    $(this).show();
                }
            });
        }
    }));
    $("#custombcsearch").keyup($.debounce(250, function() {
        const needle = $(this).val().toLowerCase();
        let hidden = 0, shown = 0;
        if (!needle) {
            $('section[id^=boundary-conditions-] section table tbody').children("tr").each(function() {
                $(this).show();
            });
            $('section[id^=boundary-conditions-] section').each(function() {
                $(this).show();
            });
            $('section[id^=boundary-conditions-]').each(function() {
                $(this).show();
            });
        } else {
            $('section[id^=boundary-conditions-] section table tbody').children("tr").each(function() {
                let found = false;
                $(this).children("td").each(function(idx, el) {
                    if ($(this).text().toLowerCase().indexOf(needle) > -1) {
                        found = true;
                    }
                });
                if (found) {
                    $(this).show();
                    $(this).parents('section[id^=boundary-conditions-] section').show();
                    $(this).parents('section[id^=boundary-conditions-]').show();
                    shown += 1;
                } else {
                    $(this).hide();
                    hidden += 1;
                }  
            });
            $('section[id^=boundary-conditions-]').each(function(idx, el) {
                $(`#${el.id} section`).each(function() {
                    if ($(this).find("table tbody tr td :visible").length === 0) {
                        $(this).hide();
                    } else {
                        $(this).show();
                    }
                });
                if ($(`#${el.id} section :visible`).length === 0) {
                    $(this).hide();
                } else {
                    $(this).show();
                }
            });
        }
    }));
    $("#custommodsearch").keyup($.debounce(250, function () {
        const needle = $(this).val().toLowerCase();
        let hidden = 0, shown = 0;
        if (!needle) {
            $('#fortran-modules section table tbody').children("tr").each(function () {
                $(this).show();
            });
            $('#fortran-modules section').each(function (idx, el) {
                $(this).show();
            });
        } else {
            $('#fortran-modules section table tbody').children("tr").each(function () {
                let found = false;
                $(this).children("td").each(function (idx, el) {
                    if ($(this).text().toLowerCase().indexOf(needle) > -1) {
                        found = true;
                    }
                });
                if (found) {
                    $(this).show();
                    $(this).parents('#fortran-modules section').show();
                    shown += 1;
                } else {
                    $(this).hide();
                    hidden += 1;
                }
            });
            $('#fortran-modules section').each(function (idx, el) {
                if ($(`#${el.id} table tbody tr td :visible`).length === 0) {
                    $(this).hide();
                } else {
                    $(this).show();
                }
            });
        }
    }));
});

/*!
 * jQuery scrollintoview() plugin and :scrollable selector filter
 *
 * Version 1.8 (14 Jul 2011)
 * Requires jQuery 1.4 or newer
 *
 * Copyright (c) 2011 Robert Koritnik
 * Licensed under the terms of the MIT license
 * http://www.opensource.org/licenses/mit-license.php
 */
(function ($) {
    var converter = {
        vertical: { x: false, y: true },
        horizontal: { x: true, y: false },
        both: { x: true, y: true },
        x: { x: true, y: false },
        y: { x: false, y: true }
    };
    var settings = {
        duration: "fast",
        direction: "both"
    };
    var rootrx = /^(?:html)$/i;
    // gets border dimensions
    var borders = function (domElement, styles) {
        styles = styles || (document.defaultView && document.defaultView.getComputedStyle ? document.defaultView.getComputedStyle(domElement, null) : domElement.currentStyle);
        var px = document.defaultView && document.defaultView.getComputedStyle ? true : false;
        var b = {
            top: (parseFloat(px ? styles.borderTopWidth : $.css(domElement, "borderTopWidth")) || 0),
            left: (parseFloat(px ? styles.borderLeftWidth : $.css(domElement, "borderLeftWidth")) || 0),
            bottom: (parseFloat(px ? styles.borderBottomWidth : $.css(domElement, "borderBottomWidth")) || 0),
            right: (parseFloat(px ? styles.borderRightWidth : $.css(domElement, "borderRightWidth")) || 0)
        };
        return {
            top: b.top,
            left: b.left,
            bottom: b.bottom,
            right: b.right,
            vertical: b.top + b.bottom,
            horizontal: b.left + b.right
        };
    };
    var dimensions = function ($element) {
        var win = $(window);
        var isRoot = rootrx.test($element[0].nodeName);
        return {
            border: isRoot ? { top: 0, left: 0, bottom: 0, right: 0 } : borders($element[0]),
            scroll: {
                top: (isRoot ? win : $element).scrollTop(),
                left: (isRoot ? win : $element).scrollLeft()
            },
            scrollbar: {
                right: isRoot ? 0 : $element.innerWidth() - $element[0].clientWidth,
                bottom: isRoot ? 0 : $element.innerHeight() - $element[0].clientHeight
            },
            rect: (function () {
                var r = $element[0].getBoundingClientRect();
                return {
                    top: isRoot ? 0 : r.top,
                    left: isRoot ? 0 : r.left,
                    bottom: isRoot ? $element[0].clientHeight : r.bottom,
                    right: isRoot ? $element[0].clientWidth : r.right
                };
            })()
        };
    };
    $.fn.extend({
        scrollintoview: function (options) {
            /// <summary>Scrolls the first element in the set into view by scrolling its closest scrollable parent.</summary>
            /// <param name="options" type="Object">Additional options that can configure scrolling:
            ///        duration (default: "fast") - jQuery animation speed (can be a duration string or number of milliseconds)
            ///        direction (default: "both") - select possible scrollings ("vertical" or "y", "horizontal" or "x", "both")
            ///        complete (default: none) - a function to call when scrolling completes (called in context of the DOM element being scrolled)
            /// </param>
            /// <return type="jQuery">Returns the same jQuery set that this function was run on.</return>
            options = $.extend({}, settings, options);
            options.direction = converter[typeof (options.direction) === "string" && options.direction.toLowerCase()] || converter.both;
            var dirStr = "";
            if (options.direction.x === true) dirStr = "horizontal";
            if (options.direction.y === true) dirStr = dirStr ? "both" : "vertical";
            var el = this.eq(0);
            var scroller = el.closest(":scrollable(" + dirStr + ")");
            // check if there's anything to scroll in the first place
            if (scroller.length > 0) {
                scroller = scroller.eq(0);
                var dim = {
                    e: dimensions(el),
                    s: dimensions(scroller)
                };
                var rel = {
                    top: dim.e.rect.top - (dim.s.rect.top + dim.s.border.top),
                    bottom: dim.s.rect.bottom - dim.s.border.bottom - dim.s.scrollbar.bottom - dim.e.rect.bottom,
                    left: dim.e.rect.left - (dim.s.rect.left + dim.s.border.left),
                    right: dim.s.rect.right - dim.s.border.right - dim.s.scrollbar.right - dim.e.rect.right
                };
                var animOptions = {};
                // vertical scroll
                if (options.direction.y === true) {
                    if (rel.top < 0) {
                        animOptions.scrollTop = dim.s.scroll.top + rel.top;
                    }
                    else if (rel.top > 0 && rel.bottom < 0) {
                        animOptions.scrollTop = dim.s.scroll.top + Math.min(rel.top, -rel.bottom);
                    }
                }
                // horizontal scroll
                if (options.direction.x === true) {
                    if (rel.left < 0) {
                        animOptions.scrollLeft = dim.s.scroll.left + rel.left;
                    }
                    else if (rel.left > 0 && rel.right < 0) {
                        animOptions.scrollLeft = dim.s.scroll.left + Math.min(rel.left, -rel.right);
                    }
                }
                // scroll if needed
                if (!$.isEmptyObject(animOptions)) {
                    if (rootrx.test(scroller[0].nodeName)) {
                        scroller = $("html,body");
                    }
                    scroller
                        .animate(animOptions, options.duration)
                        .eq(0) // we want function to be called just once (ref. "html,body")
                        .queue(function (next) {
                            $.isFunction(options.complete) && options.complete.call(scroller[0]);
                            next();
                        });
                }
                else {
                    // when there's nothing to scroll, just call the "complete" function
                    $.isFunction(options.complete) && options.complete.call(scroller[0]);
                }
            }
            // return set back
            return this;
        }
    });
    var scrollValue = {
        auto: true,
        scroll: true,
        visible: false,
        hidden: false
    };
    $.extend($.expr[":"], {
        scrollable: function (element, index, meta, stack) {
            var direction = converter[typeof (meta[3]) === "string" && meta[3].toLowerCase()] || converter.both;
            var styles = (document.defaultView && document.defaultView.getComputedStyle ? document.defaultView.getComputedStyle(element, null) : element.currentStyle);
            var overflow = {
                x: scrollValue[styles.overflowX.toLowerCase()] || false,
                y: scrollValue[styles.overflowY.toLowerCase()] || false,
                isRoot: rootrx.test(element.nodeName)
            };
            // check if completely unscrollable (exclude HTML element because it's special)
            if (!overflow.x && !overflow.y && !overflow.isRoot) {
                return false;
            }
            var size = {
                height: {
                    scroll: element.scrollHeight,
                    client: element.clientHeight
                },
                width: {
                    scroll: element.scrollWidth,
                    client: element.clientWidth
                },
                // check overflow.x/y because iPad (and possibly other tablets) don't dislay scrollbars
                scrollableX: function () {
                    return (overflow.x || overflow.isRoot) && this.width.scroll > this.width.client;
                },
                scrollableY: function () {
                    return (overflow.y || overflow.isRoot) && this.height.scroll > this.height.client;
                }
            };
            return direction.y && size.scrollableY() || direction.x && size.scrollableX();
        }
    });
})(jQuery);