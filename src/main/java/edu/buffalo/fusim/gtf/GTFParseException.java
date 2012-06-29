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
 * Exception class for GTF Parsing errors
 * 
 * @author Andrew E. Bruno
 *
 */
public class GTFParseException extends Exception {
    private static final long serialVersionUID = -7344433708122452513L;

    public GTFParseException() {
    }

    public GTFParseException(String message) {
        super(message);
    }

    public GTFParseException(Throwable cause) {
        super(cause);
    }

    public GTFParseException(String message, Throwable cause) {
        super(message, cause);
    }

}
