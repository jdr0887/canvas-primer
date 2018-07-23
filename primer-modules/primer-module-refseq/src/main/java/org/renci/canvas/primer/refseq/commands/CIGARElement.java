package org.renci.canvas.primer.refseq.commands;

import java.io.Serializable;

public class CIGARElement implements Serializable {

    private static final long serialVersionUID = -698042843410758622L;

    private CIGARType type;

    private Integer length;

    public CIGARElement() {
        super();
    }

    public CIGARElement(CIGARType type, Integer length) {
        super();
        this.type = type;
        this.length = length;
    }

    public CIGARType getType() {
        return type;
    }

    public void setType(CIGARType type) {
        this.type = type;
    }

    public Integer getLength() {
        return length;
    }

    public void setLength(Integer length) {
        this.length = length;
    }

}
