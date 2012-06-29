/*
 * Copyright 2012 Andrew E. Bruno <aebruno2@buffalo.edu>
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

package edu.buffalo.fusim.gtf;

/**
 * Enum for encapsulating a strand {+,-}
 * 
 * @author Andrew E. Bruno
 * 
 */
public enum Strand {
    FORWARD {
        public String toString() {
            return "+";
        }
    },
    REVERSE {
        public String toString() {
            return "-";
        }
    };
    
    public static Strand fromString(String str) throws GTFParseException {
        if(str.equalsIgnoreCase(Strand.FORWARD.toString())) {
            return Strand.FORWARD;
        } else if(str.equalsIgnoreCase(Strand.REVERSE.toString())) {
            return Strand.REVERSE;
        } else {
            throw new GTFParseException("Invalid strand '"+str+"'");
        }
    }
}
