<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">

<xsl:for-each select="lbe3d/core">
if (<xsl:value-of select="@name"/>)
   message(STATUS "<xsl:value-of select="@name"/> is ON")
   add_definitions(-D<xsl:value-of select="@name"/>)
   SET(COMMON ${COMMON} <xsl:value-of select="@file"/>)
   <xsl:value-of select="cmake/."/>
else(<xsl:value-of select="@name"/>)
   message(STATUS "<xsl:value-of select="@name"/> is OFF")
endif(<xsl:value-of select="@name"/>)
<!--    get_parameter_<xsl:value-of select="@type"/>(<xsl:value-of select="@name"/>,<xsl:text disable-output-escaping="yes">&amp;</xsl:text><xsl:value-of select="@name"/>);
    <xsl:for-each select="@*">
      attribute name: <xsl:value-of select="name()"/>
      attribute value: <xsl:value-of select="."/>
    </xsl:for-each>
-->
</xsl:for-each>

<!--
<xsl:for-each select="lbe3d/core">
#ifdef <xsl:value-of select="@name"/>
<xsl:for-each select="variable">
    get_parameter_<xsl:value-of select="@type"/>(<xsl:value-of select="@name"/>,<xsl:text disable-output-escaping="yes">&amp;</xsl:text><xsl:value-of select="@name"/>);
    <xsl:for-each select="@*">
      attribute name: <xsl:value-of select="name()"/>
      attribute value: <xsl:value-of select="."/>
    </xsl:for-each>
</xsl:for-each>
#endif</xsl:for-each>
-->
</xsl:template>
</xsl:stylesheet>
